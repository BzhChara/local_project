import argparse
import os
import threading
from collections import deque

import serial

from config_loader import merge_config
from excel_writer import (
    append_excel_row,
    check_output_available,
    create_raw_excel,
    create_run_output_dir,
    create_stable_excel,
    refresh_stable_excel,
    resolve_output_path,
)
from filter_logic import (
    STABLE_COUNT,
    STABLE_RANGE,
    WAVELENGTH_CODES,
    WAVELENGTH_COLORS,
    ZERO_VALUE_EPSILON,
    channel_group_label,
    channel_indices,
    collect_air_baseline,
    print_air_baseline_summary,
)
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


def parse_args():
    parser = argparse.ArgumentParser(description="持续记录分光模块串口数据")
    parser.add_argument("--port", default=None, help="串口号；未填写时读取 config.ini")
    parser.add_argument("--wavelength", type=int, choices=sorted(WAVELENGTH_CODES), default=None, help="波长；未填写时读取 config.ini")
    parser.add_argument("--channel-group", type=int, choices=(0, 1, 2), default=None, help="测试通道组：0=CH1-CH24，1=CH1-CH12，2=CH13-CH24")
    parser.add_argument("--output", default=None, help="输出 Excel 文件名；未填写时读取 config.ini")
    parser.add_argument("--stable-output", default=None, help="稳定值 Excel 文件名；未填写时读取 config.ini")
    parser.add_argument("--filter-ratio", type=int, choices=(0, 10, 20, 30), default=None, help="滤光片透过率百分比；0 表示只记录空气值")
    parser.add_argument("--ratio-tolerance", type=float, default=None, help="滤光片比例允许误差百分点；未填写时读取 config.ini")
    parser.add_argument("--air-tolerance", type=float, default=None, help="回到空气本底的允许误差百分点；未填写时读取 config.ini")
    parser.add_argument("--no-start", dest="no_start", action="store_true", default=None, help="不发送开灯命令，只读取已有数据流")
    parser.add_argument("--send-start", dest="no_start", action="store_false", default=None, help="临时覆盖 config.ini，发送开灯命令")
    parser.add_argument("--keep-light", dest="keep_light", action="store_true", default=None, help="结束时不发送关闭命令")
    parser.add_argument("--close-light", dest="keep_light", action="store_false", default=None, help="临时覆盖 config.ini，结束时发送关闭命令")
    return merge_config(parser.parse_args(), parser)


def print_startup_info(args, output_dir: str, raw_output_path: str, stable_output_path: str, target_ratio):
    if args.config_loaded:
        print(f"配置文件: {os.path.abspath(args.config_path)}")
    else:
        print(f"配置文件: 未找到 {os.path.abspath(args.config_path)}，使用命令行参数和程序默认值")

    print(f"串口: {args.port}, 波长: {args.wavelength}nm")
    print(f"串口配置: {SERIAL_BAUD} 8N1")
    print(f"预期灯色: {WAVELENGTH_COLORS[args.wavelength]}光")
    print(f"测试通道: {channel_group_label(args.channel_group)}")
    print(f"本次记录目录: {os.path.abspath(output_dir)}")
    print(f"原始数据文件: {os.path.abspath(raw_output_path)}")

    if target_ratio is not None:
        print(f"稳定值文件: {os.path.abspath(stable_output_path)}")
        print(f"稳定判定: 最近 {STABLE_COUNT} 帧极差 < {STABLE_RANGE}")
        print(f"滤光片判定: 当前值 / 空气本底 = {args.filter_ratio:g}% ± {args.ratio_tolerance:g} 个百分点")
        print(f"回到空气判定: 当前值 / 空气本底 = 100% ± {args.air_tolerance:g} 个百分点，且稳定 {STABLE_COUNT} 帧")
        print("请保持空气状态，程序将先建立目标通道空气基底")
    else:
        print("记录模式: 只记录空气值，每收到一帧有效数据就追加写入原始数据文件")


def main():
    args = parse_args()
    buffer = bytearray()
    stop_event = threading.Event()
    target_ratio = args.filter_ratio / 100 if args.filter_ratio != 0 else None
    ratio_tolerance = args.ratio_tolerance / 100
    air_tolerance = args.air_tolerance / 100
    target_indices = channel_indices(args.channel_group)
    output_dir = create_run_output_dir()
    raw_output_path = resolve_output_path(output_dir, args.output)
    stable_output_path = resolve_output_path(output_dir, args.stable_output)

    wavelength_code = WAVELENGTH_CODES[args.wavelength]
    start_cmd = build_cmd(1, 0x50, (wavelength_code << 8) | 0x01)
    stop_cmd = build_cmd(1, 0x50, 0)

    if target_ratio is not None and os.path.abspath(raw_output_path) == os.path.abspath(stable_output_path):
        print("原始数据文件和稳定值文件不能使用同一个路径")
        return

    if not check_output_available(raw_output_path):
        return
    if target_ratio is not None and not check_output_available(stable_output_path):
        return

    print_startup_info(args, output_dir, raw_output_path, stable_output_path, target_ratio)

    wb, ws = create_raw_excel(raw_output_path)
    stable_wb = stable_ws = None
    stable_windows = air_return_windows = None
    stable_values = waiting_for_air = filter_detected = None
    baseline_values = None

    if target_ratio is not None:
        stable_wb, stable_ws = create_stable_excel(stable_output_path)
        stable_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
        air_return_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
        stable_values = [None for _ in range(24)]
        waiting_for_air = [False for _ in range(24)]
        filter_detected = [False for _ in range(24)]

    row_count = 0

    with serial.Serial(
        args.port,
        SERIAL_BAUD,
        bytesize=SERIAL_BYTESIZE,
        parity=SERIAL_PARITY,
        stopbits=SERIAL_STOPBITS,
        timeout=0.1,
    ) as ser:
        ser.reset_input_buffer()

        if not args.no_start:
            ser.write(start_cmd)
            print(f"已发送启动命令: {frame_to_hex(start_cmd)}")

        if target_ratio is not None:
            baseline_values = collect_air_baseline(ser, target_indices)
            print_air_baseline_summary(baseline_values, target_indices)
            ser.reset_input_buffer()
            buffer.clear()
            print("空气基底采集完成，已清空基底阶段数据；现在可以插入滤光片")

        def wait_for_exit():
            if target_ratio is not None:
                input("开始滤光片测试，可插入或拔出滤光片；按回车退出并保存结果...")
            else:
                input("开始持续采样，按回车退出并保存结果...")
            stop_event.set()

        threading.Thread(target=wait_for_exit, daemon=True).start()

        try:
            while not stop_event.is_set():
                data = ser.read(4096)
                if data:
                    buffer.extend(data)

                for frame in parse_frames(buffer):
                    values = parse_absorbance_frame(frame)
                    if values is None:
                        continue

                    append_excel_row(wb, ws, raw_output_path, values)
                    row_count += 1

                    if target_ratio is None:
                        print(f"\r已写入 {row_count} 行", end="")
                        continue

                    stable_changed = False
                    for idx in target_indices:
                        value = values[idx]
                        if value <= ZERO_VALUE_EPSILON:
                            stable_windows[idx].clear()
                            air_return_windows[idx].clear()
                            filter_detected[idx] = False
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
                                stable_windows[idx].clear()
                                air_return_windows[idx].clear()
                                print(f"\nCH{idx + 1} 已回到空气，可再次插入滤光片")
                            continue

                        if abs(ratio - target_ratio) > ratio_tolerance:
                            stable_windows[idx].clear()
                            filter_detected[idx] = False
                            continue

                        if not filter_detected[idx]:
                            filter_detected[idx] = True
                            print(f"\nCH{idx + 1} 检测到滤光片: 当前值 {value:.6f}, 空气本底 {baseline_values[idx]:.6f}, 比例 {ratio * 100:.2f}%")

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
                            print(f"\nCH{idx + 1} 滤光片稳定: {final_value:.6f}, 比例: {final_value / baseline_values[idx] * 100:.2f}%")

                    if stable_changed:
                        refresh_stable_excel(stable_wb, stable_ws, stable_output_path, stable_values, target_indices)

                    stable_count = sum(value is not None for value in stable_values)
                    print(f"\r已写入 {row_count} 行，稳定 {stable_count}/24", end="")

        finally:
            if not args.keep_light and not args.no_start:
                ser.write(stop_cmd)
                print(f"\n已发送关闭命令: {frame_to_hex(stop_cmd)}")

    print(f"\n原始数据已保存至: {os.path.abspath(raw_output_path)}")
    if target_ratio is not None:
        refresh_stable_excel(stable_wb, stable_ws, stable_output_path, stable_values, target_indices)
        print(f"稳定值已保存至: {os.path.abspath(stable_output_path)}")


if __name__ == "__main__":
    main()
