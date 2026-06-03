import configparser
import os

from app_paths import get_app_dir


CONFIG_FILE_NAME = "config.ini"

DEFAULT_OPTIONS = {
    "port": "COM6",
    "wavelength": 520,
    "channel_group": 0,
    "output": "串口数据记录.xlsx",
    "stable_output": "滤光片稳定值.xlsx",
    "filter_ratio": None,
    "ratio_tolerance": 5.0,
    "air_tolerance": 5.0,
    "no_start": False,
    "keep_light": False,
}


def get_config_path():
    return os.path.join(get_app_dir(), CONFIG_FILE_NAME)


def load_config(path: str):
    if not os.path.exists(path):
        return {}

    with open(path, "r", encoding="utf-8") as file:
        text = file.read()

    if "[" not in text:
        text = "[settings]\n" + text

    parser = configparser.ConfigParser()
    parser.read_string(text)
    if not parser.has_section("settings"):
        return {}

    return {key.replace("-", "_"): value for key, value in parser.items("settings")}


def merge_config(args, parser):
    config_path = get_config_path()
    config_values = load_config(config_path)

    for name, default_value in DEFAULT_OPTIONS.items():
        cli_value = getattr(args, name)
        if cli_value is not None:
            value = cli_value
        elif name in config_values:
            value = convert_config_value(name, config_values[name], parser)
        else:
            value = default_value

        setattr(args, name, value)

    validate_args(args, parser)
    args.config_path = config_path
    args.config_loaded = bool(config_values)
    return args


def convert_config_value(name: str, value: str, parser):
    try:
        if name in ("wavelength", "channel_group", "filter_ratio"):
            return int(value)
        if name in ("ratio_tolerance", "air_tolerance"):
            return float(value)
        if name in ("no_start", "keep_light"):
            return parse_bool(value)
    except ValueError:
        parser.error(f"config.ini 中 {name} 的值无效: {value}")

    return value


def parse_bool(value: str):
    normalized = value.strip().lower()
    if normalized in ("1", "true", "yes", "on", "y"):
        return True
    if normalized in ("0", "false", "no", "off", "n"):
        return False
    raise ValueError(value)


def validate_args(args, parser):
    if args.wavelength not in (410, 460, 520, 550, 590, 630):
        parser.error("wavelength 必须是 410/460/520/550/590/630")
    if args.channel_group not in (0, 1, 2):
        parser.error("channel_group 必须是 0/1/2")
    if args.filter_ratio not in (0, 10, 20, 30):
        parser.error("filter_ratio 必须在 config.ini 中设置为 0/10/20/30")
