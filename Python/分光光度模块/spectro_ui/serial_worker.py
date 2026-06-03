from __future__ import annotations

import os
import shutil
from collections import deque
from dataclasses import dataclass

import serial
from openpyxl.utils.exceptions import InvalidFileException
from PySide6.QtCore import QThread, Signal
from serial.tools import list_ports

from excel_writer import (
    append_raw_csv_row,
    check_output_available,
    create_raw_csv,
    create_run_output_dir,
    create_stable_excel,
    raw_csv_to_excel,
    refresh_stable_excel,
    resolve_output_path,
)
from filter_logic import STABLE_COUNT, STABLE_RANGE, WAVELENGTH_CODES, ZERO_VALUE_EPSILON, channel_indices
from serial_protocol import (
    SERIAL_BAUD,
    SERIAL_BYTESIZE,
    SERIAL_PARITY,
    SERIAL_STOPBITS,
    build_cmd,
    frame_to_hex,
    parse_absorbance_frame,
    parse_frames,
)


def list_serial_port_names() -> list[str]:
    return [port.device for port in list_ports.comports()]


class SaveOutputError(OSError):
    """Excel 输出保存失败。"""


@dataclass(frozen=True)
class TestConfig:
    port: str
    wavelength: int
    channel_group: int
    filter_ratio: int
    output: str
    stable_output: str
    ratio_tolerance: float
    air_tolerance: float
    no_start: bool
    keep_light: bool


class SerialTestWorker(QThread):
    status_changed = Signal(str)
    log_added = Signal(str)
    values_received = Signal(object)
    baseline_progress_changed = Signal(int, int)
    baseline_value_changed = Signal(int, float)
    channel_state_changed = Signal(int, str)
    stable_value_changed = Signal(int, float, float)
    row_status_changed = Signal(int, int)
    output_paths_ready = Signal(str, str)
    test_finished = Signal()

    def __init__(self, config: TestConfig, parent=None) -> None:
        super().__init__(parent)
        self.config = config
        self._running = False

    def stop(self) -> None:
        self._running = False

    def run(self) -> None:
        self._running = True
        buffer = bytearray()
        target_indices = channel_indices(self.config.channel_group)
        target_ratio = self.config.filter_ratio / 100 if self.config.filter_ratio != 0 else None
        ratio_tolerance = self.config.ratio_tolerance / 100
        air_tolerance = self.config.air_tolerance / 100
        had_error = False
        final_status = "状态：未运行"

        output_dir = create_run_output_dir()
        raw_output_path = resolve_output_path(output_dir, self.config.output)
        stable_output_path = resolve_output_path(output_dir, self.config.stable_output)
        self.output_paths_ready.emit(os.path.abspath(raw_output_path), os.path.abspath(stable_output_path))

        if target_ratio is not None and os.path.abspath(raw_output_path) == os.path.abspath(stable_output_path):
            self.status_changed.emit("状态：原始数据文件和稳定值文件不能使用同一个路径")
            self.test_finished.emit()
            return

        if not check_output_available(raw_output_path):
            self.status_changed.emit("状态：原始数据文件被占用，请关闭后重试")
            self.test_finished.emit()
            return
        if target_ratio is not None and not check_output_available(stable_output_path):
            self.status_changed.emit("状态：稳定值文件被占用，请关闭后重试")
            self.test_finished.emit()
            return

        raw_csv_file, raw_csv_writer, raw_csv_path = create_raw_csv(raw_output_path, target_indices)
        stable_wb = stable_ws = None
        stable_values = None
        baseline_values = None

        if target_ratio is not None:
            stable_wb, stable_ws = create_stable_excel(stable_output_path, target_indices)
            stable_values = [None for _ in range(24)]

        row_count = 0
        stable_count = 0
        start_cmd = build_cmd(1, 0x50, (WAVELENGTH_CODES[self.config.wavelength] << 8) | 0x01)
        stop_cmd = build_cmd(1, 0x50, 0)

        ser = None
        try:
            ser = serial.Serial(
                self.config.port,
                SERIAL_BAUD,
                bytesize=SERIAL_BYTESIZE,
                parity=SERIAL_PARITY,
                stopbits=SERIAL_STOPBITS,
                timeout=0.1,
                write_timeout=1,
            )
            ser.reset_input_buffer()
            self.status_changed.emit(f"状态：已连接 {self.config.port}")
            self.log_added.emit(f"串口配置: {SERIAL_BAUD} 8N1")

            if not self.config.no_start:
                ser.write(start_cmd)
                self.log_added.emit(f"已发送启动命令: {frame_to_hex(start_cmd)}")

            if target_ratio is not None:
                self.status_changed.emit("状态：请保持空气状态，正在建立空气基底")
                baseline_values = self._collect_air_baseline(ser, target_indices)
                if not self._running:
                    return
                ser.reset_input_buffer()
                buffer.clear()
                self.log_added.emit("空气基底采集完成，已清空基底阶段数据；现在可以插入滤光片")
                self.status_changed.emit("状态：基底完成，可以插入滤光片")
            else:
                for idx in target_indices:
                    self.channel_state_changed.emit(idx, "ready")
                self.status_changed.emit("状态：空气记录模式，正在持续采样")

            stable_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
            air_return_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
            waiting_for_air = [False for _ in range(24)]
            filter_detected = [False for _ in range(24)]

            while self._running:
                data = ser.read(4096)
                if data:
                    buffer.extend(data)

                for frame in parse_frames(buffer):
                    values = parse_absorbance_frame(frame)
                    if values is None:
                        continue

                    self.values_received.emit(values)
                    try:
                        append_raw_csv_row(raw_csv_file, raw_csv_writer, values, target_indices)
                    except (OSError, PermissionError, ValueError) as exc:
                        raise SaveOutputError(f"原始 CSV 缓存保存失败: {exc}") from exc
                    row_count += 1

                    if target_ratio is None:
                        self.row_status_changed.emit(row_count, stable_count)
                        continue

                    stable_changed = False
                    for idx in target_indices:
                        value = values[idx]
                        if value <= ZERO_VALUE_EPSILON:
                            stable_windows[idx].clear()
                            air_return_windows[idx].clear()
                            if waiting_for_air[idx]:
                                continue
                            if filter_detected[idx]:
                                filter_detected[idx] = False
                                self.channel_state_changed.emit(idx, "ready")
                            continue

                        ratio = value / baseline_values[idx]
                        if waiting_for_air[idx]:
                            if abs(ratio - 1.0) > air_tolerance:
                                air_return_windows[idx].clear()
                                continue

                            air_return_windows[idx].append(value)
                            if len(air_return_windows[idx]) < STABLE_COUNT:
                                continue

                            air_range = max(air_return_windows[idx]) - min(air_return_windows[idx])
                            if air_range < STABLE_RANGE:
                                waiting_for_air[idx] = False
                                filter_detected[idx] = False
                                stable_windows[idx].clear()
                                air_return_windows[idx].clear()
                                self.channel_state_changed.emit(idx, "ready")
                                self.log_added.emit(f"CH{idx + 1} 已回到空气，可再次插入滤光片")
                            continue

                        if abs(ratio - target_ratio) > ratio_tolerance:
                            stable_windows[idx].clear()
                            if filter_detected[idx]:
                                filter_detected[idx] = False
                                self.channel_state_changed.emit(idx, "ready")
                            continue

                        if not filter_detected[idx]:
                            filter_detected[idx] = True
                            self.channel_state_changed.emit(idx, "detecting")
                            self.log_added.emit(
                                f"CH{idx + 1} 检测到滤光片: 当前值 {value:.6f}, "
                                f"空气本底 {baseline_values[idx]:.6f}, 比例 {ratio * 100:.2f}%"
                            )

                        stable_windows[idx].append(value)
                        if len(stable_windows[idx]) < STABLE_COUNT:
                            continue

                        sample_range = max(stable_windows[idx]) - min(stable_windows[idx])
                        if sample_range < STABLE_RANGE:
                            final_value = round(sum(stable_windows[idx]) / STABLE_COUNT, 6)
                            stable_values[idx] = final_value
                            waiting_for_air[idx] = True
                            stable_windows[idx].clear()
                            air_return_windows[idx].clear()
                            stable_changed = True
                            ratio_percent = final_value / baseline_values[idx] * 100
                            self.channel_state_changed.emit(idx, "waiting_air")
                            self.stable_value_changed.emit(idx, final_value, ratio_percent)
                            self.log_added.emit(
                                f"CH{idx + 1} 滤光片稳定: {final_value:.6f}, 比例: {ratio_percent:.2f}%；请拔出滤光片"
                            )

                    if stable_changed:
                        stable_count = sum(value is not None for value in stable_values)
                        try:
                            refresh_stable_excel(stable_wb, stable_ws, stable_output_path, stable_values, target_indices)
                        except (OSError, PermissionError, InvalidFileException) as exc:
                            raise SaveOutputError(f"稳定值保存失败: {exc}") from exc

                    self.row_status_changed.emit(row_count, stable_count)

        except SaveOutputError as exc:
            had_error = True
            final_status = f"状态：保存异常：{exc}"
            if self._running:
                self.status_changed.emit(final_status)
                self.log_added.emit(f"保存异常: {exc}")
                self.log_added.emit("请确认 Excel/WPS 没有打开输出文件，资源管理器未预览该文件，并避免杀毒或同步软件占用 data 目录")
        except Exception as exc:
            had_error = True
            final_status = f"状态：串口异常：{exc}"
            if self._running:
                self.status_changed.emit(final_status)
                self.log_added.emit(f"串口异常: {exc}")
        finally:
            self._running = False
            try:
                if ser is not None and ser.is_open and not self.config.keep_light and not self.config.no_start:
                    ser.write(stop_cmd)
                    self.log_added.emit(f"已发送关闭命令: {frame_to_hex(stop_cmd)}")
            except Exception as exc:
                self.log_added.emit(f"发送关闭命令失败: {exc}")
            finally:
                if ser is not None and ser.is_open:
                    ser.close()

            try:
                raw_csv_file.close()
            except OSError as exc:
                self.log_added.emit(f"关闭原始 CSV 缓存失败: {exc}")

            should_keep_raw = row_count > 0
            should_keep_stable = target_ratio is not None and stable_count > 0
            raw_excel_saved = False

            if should_keep_raw:
                try:
                    raw_csv_to_excel(raw_csv_path, raw_output_path)
                    raw_excel_saved = True
                except (OSError, PermissionError, ValueError, InvalidFileException) as exc:
                    had_error = True
                    final_status = f"状态：保存异常：原始 Excel 生成失败: {exc}"
                    self.log_added.emit(f"原始 Excel 生成失败: {exc}")

            if target_ratio is not None and stable_wb is not None and should_keep_stable:
                try:
                    refresh_stable_excel(stable_wb, stable_ws, stable_output_path, stable_values, target_indices)
                except (OSError, PermissionError, InvalidFileException) as exc:
                    self.log_added.emit(f"保存稳定值失败: {exc}")

            if should_keep_raw and raw_excel_saved:
                self.log_added.emit(f"原始数据已保存至: {os.path.abspath(raw_output_path)}")
            if should_keep_raw:
                self.log_added.emit(f"原始 CSV 缓存已保留: {os.path.abspath(raw_csv_path)}")
            if target_ratio is not None and should_keep_stable:
                self.log_added.emit(f"稳定值已保存至: {os.path.abspath(stable_output_path)}")

            if had_error or not should_keep_raw or (target_ratio is not None and not should_keep_stable):
                self._cleanup_empty_outputs(
                    output_dir,
                    raw_output_path,
                    stable_output_path if target_ratio is not None else None,
                    keep_raw=should_keep_raw,
                    keep_stable=should_keep_stable,
                    raw_csv_path=raw_csv_path,
                )
            self.status_changed.emit(final_status)
            self.test_finished.emit()

    def _cleanup_empty_outputs(
        self,
        output_dir: str,
        raw_output_path: str,
        stable_output_path: str | None,
        keep_raw: bool,
        keep_stable: bool,
        raw_csv_path: str | None = None,
    ) -> None:
        removed_files = []

        if not keep_raw and self._remove_file(raw_output_path):
            removed_files.append(os.path.basename(raw_output_path))
        if raw_csv_path is not None and not keep_raw and self._remove_file(raw_csv_path):
            removed_files.append(os.path.basename(raw_csv_path))

        if stable_output_path is not None and not keep_stable and self._remove_file(stable_output_path):
            removed_files.append(os.path.basename(stable_output_path))

        if removed_files:
            self.log_added.emit(f"未产生有效数据，已删除空文件: {', '.join(removed_files)}")

        try:
            temp_csv_dir = os.path.join(output_dir, "_temp_csv")
            if os.path.isdir(temp_csv_dir) and not os.listdir(temp_csv_dir):
                shutil.rmtree(temp_csv_dir)
            if os.path.isdir(output_dir) and not os.listdir(output_dir):
                shutil.rmtree(output_dir)
                self.log_added.emit(f"本次记录目录为空，已删除: {os.path.abspath(output_dir)}")
        except OSError as exc:
            self.log_added.emit(f"清理空记录目录失败: {exc}")

    def _remove_file(self, path: str) -> bool:
        try:
            if os.path.exists(path):
                os.remove(path)
                return True
        except OSError as exc:
            self.log_added.emit(f"删除空文件失败: {os.path.abspath(path)}: {exc}")
        return False

    def _collect_air_baseline(self, ser, target_indices: list[int]):
        buffer = bytearray()
        baseline_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
        baseline_values = [None for _ in range(24)]
        last_done_count = -1

        self.baseline_progress_changed.emit(0, len(target_indices))
        while self._running:
            data = ser.read(4096)
            if data:
                buffer.extend(data)

            for frame in parse_frames(buffer):
                values = parse_absorbance_frame(frame)
                if values is None:
                    continue

                self.values_received.emit(values)
                for idx in target_indices:
                    if baseline_values[idx] is not None:
                        continue

                    value = values[idx]
                    if value <= ZERO_VALUE_EPSILON:
                        baseline_windows[idx].clear()
                        continue

                    baseline_windows[idx].append(value)
                    if len(baseline_windows[idx]) < STABLE_COUNT:
                        continue

                    baseline_range = max(baseline_windows[idx]) - min(baseline_windows[idx])
                    if baseline_range < STABLE_RANGE:
                        baseline_values[idx] = round(sum(baseline_windows[idx]) / STABLE_COUNT, 6)
                        self.baseline_value_changed.emit(idx, baseline_values[idx])
                        self.channel_state_changed.emit(idx, "ready")

                done_count = sum(baseline_values[idx] is not None for idx in target_indices)
                if done_count != last_done_count:
                    self.baseline_progress_changed.emit(done_count, len(target_indices))
                    last_done_count = done_count

                if done_count == len(target_indices):
                    return baseline_values

        return baseline_values
