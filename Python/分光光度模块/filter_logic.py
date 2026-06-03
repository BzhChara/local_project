from collections import deque

from serial_protocol import parse_absorbance_frame, parse_frames


WAVELENGTH_CODES = {
    410: 0,
    460: 1,
    520: 2,
    550: 3,
    590: 4,
    630: 5,
}

WAVELENGTH_COLORS = {
    410: "紫色",
    460: "蓝色",
    520: "绿色",
    550: "绿色",
    590: "橙色",
    630: "红色",
}

ZERO_VALUE_EPSILON = 1e-7
STABLE_COUNT = 10
STABLE_RANGE = 0.005


def channel_indices(channel_group: int):
    if channel_group == 1:
        return list(range(12))
    if channel_group == 2:
        return list(range(12, 24))
    return list(range(24))


def channel_group_label(channel_group: int):
    if channel_group == 1:
        return "CH1-CH12"
    if channel_group == 2:
        return "CH13-CH24"
    return "CH1-CH24"


def collect_air_baseline(ser, target_indices: list[int]):
    buffer = bytearray()
    baseline_windows = [deque(maxlen=STABLE_COUNT) for _ in range(24)]
    baseline_values = [None for _ in range(24)]
    last_done_count = -1

    print(f"空气基底采集中: 0/{len(target_indices)}", end="")
    while True:
        data = ser.read(4096)
        if data:
            buffer.extend(data)

        for frame in parse_frames(buffer):
            values = parse_absorbance_frame(frame)
            if values is None:
                continue

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

            done_count = sum(baseline_values[idx] is not None for idx in target_indices)
            if done_count != last_done_count:
                print(f"\r空气基底采集中: {done_count}/{len(target_indices)}", end="")
                last_done_count = done_count

            if done_count == len(target_indices):
                print()
                return baseline_values


def print_air_baseline_summary(baseline_values: list[float | None], target_indices: list[int]):
    print("空气基底完成:")
    for idx in target_indices:
        value = baseline_values[idx]
        if value is not None:
            print(f"CH{idx + 1}: {value:.6f}")
