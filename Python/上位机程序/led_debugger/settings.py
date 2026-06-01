from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any


SETTINGS_PATH = Path(__file__).resolve().parent.parent / "settings.json"

DEFAULT_BAUD_RATE = 115200
DEFAULT_READ_TIMEOUT_SECONDS = 0.25
DEFAULT_ROUND_DELAY_SECONDS = 0.05
DEFAULT_CHANNEL_DELAY_SECONDS = 0.5
DEFAULT_CURVE_VISIBLE_SECONDS = 600.0
DEFAULT_MAX_POINTS = 500
DEFAULT_DISPLAY_DECIMALS = 3
DEFAULT_CURVE_SHOW_GRID = True
DEFAULT_CURVE_LINE_WIDTH = 2
DEFAULT_Y_AXIS_AUTO_SCALE = True
DEFAULT_Y_AXIS_MIN = 0.0
DEFAULT_Y_AXIS_MAX = 1.0
DEFAULT_REMEMBER_LAST_PORT = False
DEFAULT_AUTO_CONNECT_LAST_PORT = False
DEFAULT_SHOW_SERIAL_ERRORS = True
DEFAULT_SHOW_RAW_SERIAL_VIEW = False
DEFAULT_SAVE_DURATION_SECONDS = 600.0
DEFAULT_SAVE_INTERVAL_SECONDS = 10.0
DEFAULT_AUTO_SAVE_INTERVAL = True
DEFAULT_AUTO_SAVE_DURATION = True
DEFAULT_AUTO_CURVE_VISIBLE_SECONDS = True
DEFAULT_FULL_SAVE_MODE = False

READ_TIMEOUT_RANGE = (0.05, 2.0)
ROUND_DELAY_RANGE = (0.0, 1.0)
CHANNEL_DELAY_RANGE = (0.0, 5.0)
CURVE_VISIBLE_SECONDS_RANGE = (10.0, 3600.0)
MAX_POINTS_RANGE = (100, 10000)
DISPLAY_DECIMALS_RANGE = (0, 6)
CURVE_LINE_WIDTH_RANGE = (1.0, 6.0)
Y_AXIS_RANGE = (-1000000.0, 1000000.0)
SAVE_DURATION_RANGE = (1.0, 86400.0)
SAVE_INTERVAL_RANGE = (0.1, 3600.0)
SAVE_INTERVAL_CHANNEL_COUNT = 16
AUTO_DURATION_TARGET_POINTS = 100


@dataclass(frozen=True)
class AppSettings:
    """用户可调整的应用参数。"""

    read_timeout_seconds: float = DEFAULT_READ_TIMEOUT_SECONDS
    round_delay_seconds: float = DEFAULT_ROUND_DELAY_SECONDS
    channel_delay_seconds: float = DEFAULT_CHANNEL_DELAY_SECONDS
    curve_visible_seconds: float = DEFAULT_CURVE_VISIBLE_SECONDS
    max_points: int = DEFAULT_MAX_POINTS
    remember_last_port: bool = DEFAULT_REMEMBER_LAST_PORT
    last_port_name: str = ""
    auto_connect_last_port: bool = DEFAULT_AUTO_CONNECT_LAST_PORT
    display_decimals: int = DEFAULT_DISPLAY_DECIMALS
    curve_show_grid: bool = DEFAULT_CURVE_SHOW_GRID
    curve_line_width: float = DEFAULT_CURVE_LINE_WIDTH
    y_axis_auto_scale: bool = DEFAULT_Y_AXIS_AUTO_SCALE
    y_axis_min: float = DEFAULT_Y_AXIS_MIN
    y_axis_max: float = DEFAULT_Y_AXIS_MAX
    show_serial_errors: bool = DEFAULT_SHOW_SERIAL_ERRORS
    show_raw_serial_view: bool = DEFAULT_SHOW_RAW_SERIAL_VIEW
    save_duration_seconds: float = DEFAULT_SAVE_DURATION_SECONDS
    save_interval_seconds: float = DEFAULT_SAVE_INTERVAL_SECONDS
    auto_save_interval: bool = DEFAULT_AUTO_SAVE_INTERVAL
    auto_save_duration: bool = DEFAULT_AUTO_SAVE_DURATION
    auto_curve_visible_seconds: bool = DEFAULT_AUTO_CURVE_VISIBLE_SECONDS
    full_save_mode: bool = DEFAULT_FULL_SAVE_MODE


def load_settings(path: Path = SETTINGS_PATH) -> AppSettings:
    """从本地 JSON 读取设置；文件不存在或格式异常时使用默认值。"""
    if not path.exists():
        return AppSettings()

    try:
        raw_data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return AppSettings()

    if not isinstance(raw_data, dict):
        return AppSettings()

    allowed_names = {field.name for field in fields(AppSettings)}
    values: dict[str, Any] = {}
    for name in allowed_names:
        if name in raw_data:
            values[name] = raw_data[name]

    try:
        settings = AppSettings(
            read_timeout_seconds=float(values.get("read_timeout_seconds", DEFAULT_READ_TIMEOUT_SECONDS)),
            round_delay_seconds=float(values.get("round_delay_seconds", DEFAULT_ROUND_DELAY_SECONDS)),
            channel_delay_seconds=float(values.get("channel_delay_seconds", DEFAULT_CHANNEL_DELAY_SECONDS)),
            curve_visible_seconds=float(values.get("curve_visible_seconds", DEFAULT_CURVE_VISIBLE_SECONDS)),
            max_points=int(values.get("max_points", DEFAULT_MAX_POINTS)),
            remember_last_port=bool(values.get("remember_last_port", DEFAULT_REMEMBER_LAST_PORT)),
            last_port_name=str(values.get("last_port_name", "")),
            auto_connect_last_port=bool(values.get("auto_connect_last_port", DEFAULT_AUTO_CONNECT_LAST_PORT)),
            display_decimals=int(values.get("display_decimals", DEFAULT_DISPLAY_DECIMALS)),
            curve_show_grid=bool(values.get("curve_show_grid", DEFAULT_CURVE_SHOW_GRID)),
            curve_line_width=float(values.get("curve_line_width", DEFAULT_CURVE_LINE_WIDTH)),
            y_axis_auto_scale=bool(values.get("y_axis_auto_scale", DEFAULT_Y_AXIS_AUTO_SCALE)),
            y_axis_min=float(values.get("y_axis_min", DEFAULT_Y_AXIS_MIN)),
            y_axis_max=float(values.get("y_axis_max", DEFAULT_Y_AXIS_MAX)),
            show_serial_errors=bool(values.get("show_serial_errors", DEFAULT_SHOW_SERIAL_ERRORS)),
            show_raw_serial_view=bool(values.get("show_raw_serial_view", DEFAULT_SHOW_RAW_SERIAL_VIEW)),
            save_duration_seconds=float(values.get("save_duration_seconds", DEFAULT_SAVE_DURATION_SECONDS)),
            save_interval_seconds=float(values.get("save_interval_seconds", DEFAULT_SAVE_INTERVAL_SECONDS)),
            auto_save_interval=bool(values.get("auto_save_interval", DEFAULT_AUTO_SAVE_INTERVAL)),
            auto_save_duration=bool(values.get("auto_save_duration", DEFAULT_AUTO_SAVE_DURATION)),
            auto_curve_visible_seconds=bool(
                values.get("auto_curve_visible_seconds", DEFAULT_AUTO_CURVE_VISIBLE_SECONDS)
            ),
            full_save_mode=bool(values.get("full_save_mode", DEFAULT_FULL_SAVE_MODE)),
        )
    except (TypeError, ValueError):
        return AppSettings()
    return normalize_settings(settings)


def save_settings(settings: AppSettings, path: Path = SETTINGS_PATH) -> None:
    """保存设置到本地 JSON 文件。"""
    path.write_text(
        json.dumps(asdict(normalize_settings(settings)), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def normalize_settings(settings: AppSettings) -> AppSettings:
    """把外部 JSON 或界面输入限制在当前软件支持范围内。"""
    return AppSettings(
        read_timeout_seconds=_clamp_float(settings.read_timeout_seconds, READ_TIMEOUT_RANGE),
        round_delay_seconds=_clamp_float(settings.round_delay_seconds, ROUND_DELAY_RANGE),
        channel_delay_seconds=_clamp_float(settings.channel_delay_seconds, CHANNEL_DELAY_RANGE),
        curve_visible_seconds=_clamp_float(settings.curve_visible_seconds, CURVE_VISIBLE_SECONDS_RANGE),
        max_points=_clamp_int(settings.max_points, MAX_POINTS_RANGE),
        remember_last_port=bool(settings.remember_last_port),
        last_port_name=str(settings.last_port_name),
        auto_connect_last_port=bool(settings.auto_connect_last_port),
        display_decimals=_clamp_int(settings.display_decimals, DISPLAY_DECIMALS_RANGE),
        curve_show_grid=bool(settings.curve_show_grid),
        curve_line_width=_clamp_float(settings.curve_line_width, CURVE_LINE_WIDTH_RANGE),
        y_axis_auto_scale=bool(settings.y_axis_auto_scale),
        y_axis_min=_clamp_float(min(settings.y_axis_min, settings.y_axis_max), Y_AXIS_RANGE),
        y_axis_max=_clamp_float(max(settings.y_axis_min, settings.y_axis_max), Y_AXIS_RANGE),
        show_serial_errors=bool(settings.show_serial_errors),
        show_raw_serial_view=bool(settings.show_raw_serial_view),
        save_duration_seconds=_clamp_float(settings.save_duration_seconds, SAVE_DURATION_RANGE),
        save_interval_seconds=_clamp_float(settings.save_interval_seconds, SAVE_INTERVAL_RANGE),
        auto_save_interval=bool(settings.auto_save_interval),
        auto_save_duration=bool(settings.auto_save_duration),
        auto_curve_visible_seconds=bool(settings.auto_curve_visible_seconds),
        full_save_mode=bool(settings.full_save_mode),
    )


def calculated_auto_save_interval_seconds(settings: AppSettings) -> float:
    """按 16 通道轮询设置估算自动保存间隔。"""
    interval_seconds = settings.channel_delay_seconds * SAVE_INTERVAL_CHANNEL_COUNT + settings.round_delay_seconds
    return _clamp_float(interval_seconds, SAVE_INTERVAL_RANGE)


def effective_save_interval_seconds(settings: AppSettings) -> float:
    """返回保存时真正使用的时间间隔。"""
    if settings.auto_save_interval:
        return calculated_auto_save_interval_seconds(settings)
    return settings.save_interval_seconds


def calculated_auto_save_duration_seconds(settings: AppSettings) -> float:
    """按实际保存间隔估算自动保存时长。"""
    return _clamp_float(effective_save_interval_seconds(settings) * AUTO_DURATION_TARGET_POINTS, SAVE_DURATION_RANGE)


def effective_save_duration_seconds(settings: AppSettings) -> float:
    """返回保存时真正使用的保存时长。"""
    if settings.auto_save_duration:
        return calculated_auto_save_duration_seconds(settings)
    return settings.save_duration_seconds


def calculated_auto_curve_visible_seconds(settings: AppSettings) -> float:
    """按实际保存间隔估算自动曲线显示时长。"""
    return _clamp_float(
        effective_save_interval_seconds(settings) * AUTO_DURATION_TARGET_POINTS,
        CURVE_VISIBLE_SECONDS_RANGE,
    )


def effective_curve_visible_seconds(settings: AppSettings) -> float:
    """返回曲线真正使用的显示时长。"""
    if settings.auto_curve_visible_seconds:
        return calculated_auto_curve_visible_seconds(settings)
    return settings.curve_visible_seconds


def _clamp_float(value: float, value_range: tuple[float, float]) -> float:
    return min(max(float(value), value_range[0]), value_range[1])


def _clamp_int(value: int, value_range: tuple[int, int]) -> int:
    return min(max(int(value), value_range[0]), value_range[1])
