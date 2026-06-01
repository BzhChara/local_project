from dataclasses import replace

import pytest

from led_debugger.settings import (
    AppSettings,
    calculated_auto_curve_visible_seconds,
    calculated_auto_save_duration_seconds,
    calculated_auto_save_interval_seconds,
    effective_curve_visible_seconds,
    effective_save_duration_seconds,
    effective_save_interval_seconds,
    load_settings,
    save_settings,
)


def test_load_settings_returns_defaults_when_file_missing(tmp_path) -> None:
    settings = load_settings(tmp_path / "missing.json")

    assert settings == AppSettings()
    assert settings.curve_visible_seconds == 600.0
    assert settings.save_duration_seconds == 600.0
    assert settings.save_interval_seconds == 10.0
    assert settings.auto_save_interval is True
    assert settings.auto_save_duration is True
    assert settings.auto_curve_visible_seconds is True
    assert settings.full_save_mode is False


def test_save_and_load_settings_roundtrip(tmp_path) -> None:
    settings_path = tmp_path / "settings.json"
    settings = AppSettings(
        read_timeout_seconds=0.5,
        round_delay_seconds=0.1,
        channel_delay_seconds=0.5,
        curve_visible_seconds=300.0,
        max_points=1000,
        auto_save_interval=False,
        auto_save_duration=False,
        auto_curve_visible_seconds=False,
        full_save_mode=True,
    )

    save_settings(settings, settings_path)

    assert load_settings(settings_path) == settings


def test_load_settings_clamps_out_of_range_values(tmp_path) -> None:
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        """
        {
          "read_timeout_seconds": 10,
          "round_delay_seconds": -1,
          "channel_delay_seconds": 10,
          "curve_visible_seconds": 1,
          "max_points": 1,
          "show_raw_serial_view": true,
          "save_duration_seconds": 0,
          "save_interval_seconds": 10000,
          "auto_save_interval": false,
          "auto_save_duration": false,
          "auto_curve_visible_seconds": false,
          "full_save_mode": true
        }
        """,
        encoding="utf-8",
    )

    settings = load_settings(settings_path)

    assert settings.read_timeout_seconds == 2.0
    assert settings.round_delay_seconds == 0.0
    assert settings.channel_delay_seconds == 5.0
    assert settings.curve_visible_seconds == 10.0
    assert settings.max_points == 100
    assert settings.show_raw_serial_view is True
    assert settings.save_duration_seconds == 1.0
    assert settings.save_interval_seconds == 3600.0
    assert settings.auto_save_interval is False
    assert settings.auto_save_duration is False
    assert settings.auto_curve_visible_seconds is False
    assert settings.full_save_mode is True


def test_effective_save_interval_uses_auto_calculation_when_enabled() -> None:
    settings = AppSettings(
        channel_delay_seconds=0.5,
        round_delay_seconds=0.05,
        save_interval_seconds=10.0,
        auto_save_interval=True,
    )

    assert calculated_auto_save_interval_seconds(settings) == 8.05
    assert effective_save_interval_seconds(settings) == 8.05
    assert effective_save_interval_seconds(replace(settings, auto_save_interval=False)) == 10.0


def test_effective_durations_use_100_point_auto_calculation_when_enabled() -> None:
    settings = AppSettings(
        channel_delay_seconds=0.5,
        round_delay_seconds=0.05,
        save_interval_seconds=10.0,
        save_duration_seconds=600.0,
        curve_visible_seconds=600.0,
        auto_save_interval=True,
        auto_save_duration=True,
        auto_curve_visible_seconds=True,
    )

    assert calculated_auto_save_duration_seconds(settings) == pytest.approx(805.0)
    assert effective_save_duration_seconds(settings) == pytest.approx(805.0)
    assert calculated_auto_curve_visible_seconds(settings) == pytest.approx(805.0)
    assert effective_curve_visible_seconds(settings) == pytest.approx(805.0)
    manual_settings = replace(settings, auto_save_duration=False, auto_curve_visible_seconds=False)
    assert effective_save_duration_seconds(manual_settings) == 600.0
    assert effective_curve_visible_seconds(manual_settings) == 600.0
