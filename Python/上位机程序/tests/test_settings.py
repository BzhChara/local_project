from led_debugger.settings import AppSettings, load_settings, save_settings


def test_load_settings_returns_defaults_when_file_missing(tmp_path) -> None:
    settings = load_settings(tmp_path / "missing.json")

    assert settings == AppSettings()


def test_save_and_load_settings_roundtrip(tmp_path) -> None:
    settings_path = tmp_path / "settings.json"
    settings = AppSettings(
        read_timeout_seconds=0.5,
        round_delay_seconds=0.1,
        curve_visible_seconds=300.0,
        max_points=1000,
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
          "curve_visible_seconds": 1,
          "max_points": 1
        }
        """,
        encoding="utf-8",
    )

    settings = load_settings(settings_path)

    assert settings.read_timeout_seconds == 2.0
    assert settings.round_delay_seconds == 0.0
    assert settings.curve_visible_seconds == 10.0
    assert settings.max_points == 100
