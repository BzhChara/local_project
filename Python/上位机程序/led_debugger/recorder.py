from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TextIO

from led_debugger.data_store import HistoryPoint
from led_debugger.exporter import (
    EXPORT_CHANNELS,
    EXPORT_LEDS,
    ExportSummary,
    application_root,
    format_export_date,
    write_history_xlsx,
)


TEMP_CSV_DIR_NAME = "_temp_csv"


@dataclass
class _TempCsvWriter:
    file: TextIO
    writer: csv.writer


def default_recording_dir(moment: datetime | None = None, data_root: Path | None = None) -> Path:
    """按开始记录时间生成连续记录目录。"""
    start_time = moment or datetime.now()
    root = data_root or application_root() / "data"
    day_dir = root / format_export_date(start_time)
    return day_dir / f"record_{start_time:%H%M%S}"


def prepare_recording_dir(moment: datetime | None = None, data_root: Path | None = None) -> Path:
    """准备连续记录目录；同一秒重复启动时覆盖旧记录。"""
    output_dir = default_recording_dir(moment, data_root)
    if not output_dir.exists():
        return output_dir
    if not output_dir.is_dir():
        raise OSError(f"记录路径已存在且不是文件夹：{output_dir}")
    _ensure_auto_recording_dir(output_dir)
    shutil.rmtree(output_dir)
    return output_dir


def _ensure_auto_recording_dir(path: Path) -> None:
    """只允许清空程序自动生成的 data/YYYY_MM_DD/record_HHMMSS 目录。"""
    name = path.name
    if len(name) != len("record_HHMMSS") or not name.startswith("record_") or not name[7:].isdigit():
        raise OSError(f"拒绝覆盖非记录目录：{path}")
    try:
        datetime.strptime(path.parent.name, "%Y_%m_%d")
    except ValueError as exc:
        raise OSError(f"拒绝覆盖日期目录异常的记录目录：{path}") from exc


class ContinuousRecorder:
    """把长时间测试数据追加到临时 CSV，停止后再生成 Excel。"""

    def __init__(self, output_dir: Path, interval_seconds: float) -> None:
        self.output_dir = output_dir
        self.temp_dir = output_dir / TEMP_CSV_DIR_NAME
        self.interval_seconds = interval_seconds
        self.rows_written = 0
        self.recorded_channel_samples = 0
        self._next_timestamp_by_channel: dict[int, float] = {}
        self._writers: dict[tuple[int, int], _TempCsvWriter] = {}
        self._closed = False
        self._open_temp_files()

    def record_channel_reading(self, channel: int, values_ma: tuple[float, ...], timestamp: float) -> int:
        """记录达到保存间隔的一组通道数据，返回新增行数。"""
        if self._closed:
            raise RuntimeError("Recorder is closed")
        if channel not in EXPORT_CHANNELS:
            return 0

        next_timestamp = self._next_timestamp_by_channel.get(channel)
        if next_timestamp is not None and timestamp < next_timestamp:
            return 0

        self._next_timestamp_by_channel[channel] = timestamp + self.interval_seconds
        rows_added = 0
        for led_index, value_ma in enumerate(values_ma[: len(EXPORT_LEDS)], start=1):
            writer = self._writers[(channel, led_index)]
            writer.writer.writerow([repr(timestamp), repr(value_ma)])
            rows_added += 1
        for led_index in EXPORT_LEDS:
            self._writers[(channel, led_index)].file.flush()

        self.rows_written += rows_added
        self.recorded_channel_samples += 1
        return rows_added

    def finish(self) -> ExportSummary:
        """关闭临时 CSV，并把全部记录转换为 Excel。"""
        self.close()
        files_written = 0
        rows_written = 0
        for channel in EXPORT_CHANNELS:
            channel_dir = self.output_dir / f"通道{channel}"
            channel_dir.mkdir(parents=True, exist_ok=True)
            for led_index in EXPORT_LEDS:
                points = self._read_temp_points(channel, led_index)
                write_history_xlsx(channel_dir / f"LED电流{led_index}.xlsx", points)
                files_written += 1
                rows_written += len(points)
        return ExportSummary(output_dir=self.output_dir, files_written=files_written, rows_written=rows_written)

    def close(self) -> None:
        """关闭临时 CSV 句柄，保留文件用于排查或恢复。"""
        if self._closed:
            return
        for temp_writer in self._writers.values():
            temp_writer.file.close()
        self._closed = True

    def _open_temp_files(self) -> None:
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            for channel in EXPORT_CHANNELS:
                for led_index in EXPORT_LEDS:
                    file = self._temp_csv_path(channel, led_index).open("w", newline="", encoding="utf-8")
                    writer = csv.writer(file)
                    writer.writerow(["timestamp", "value_ma"])
                    self._writers[(channel, led_index)] = _TempCsvWriter(file=file, writer=writer)
        except OSError:
            self.close()
            raise

    def _read_temp_points(self, channel: int, led_index: int) -> list[HistoryPoint]:
        path = self._temp_csv_path(channel, led_index)
        points: list[HistoryPoint] = []
        with path.open("r", newline="", encoding="utf-8") as file:
            reader = csv.DictReader(file)
            for row in reader:
                points.append(HistoryPoint(float(row["timestamp"]), float(row["value_ma"])))
        return points

    def _temp_csv_path(self, channel: int, led_index: int) -> Path:
        return self.temp_dir / f"通道{channel}_LED电流{led_index}.csv"
