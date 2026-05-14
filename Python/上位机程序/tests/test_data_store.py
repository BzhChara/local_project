from led_debugger.data_store import DataStore


def test_add_channel_reading_updates_latest_values() -> None:
    store = DataStore()
    values = tuple(float(index) for index in range(1, 17))

    store.add_channel_reading(channel=1, values_ma=values, timestamp=100.0)

    assert store.get_latest(1) == values


def test_add_channel_reading_appends_history_for_each_led() -> None:
    store = DataStore()
    values = tuple(float(index) for index in range(1, 17))

    store.add_channel_reading(channel=2, values_ma=values, timestamp=123.5)

    led_1_history = store.get_history(channel=2, led_index=1)
    led_16_history = store.get_history(channel=2, led_index=16)
    assert led_1_history[0].timestamp == 123.5
    assert led_1_history[0].value_ma == 1.0
    assert led_16_history[0].timestamp == 123.5
    assert led_16_history[0].value_ma == 16.0


def test_history_keeps_recent_points_only() -> None:
    store = DataStore(max_points=2)
    first = tuple(float(index) for index in range(1, 17))
    second = tuple(float(index + 10) for index in range(1, 17))
    third = tuple(float(index + 20) for index in range(1, 17))

    store.add_channel_reading(channel=1, values_ma=first, timestamp=1.0)
    store.add_channel_reading(channel=1, values_ma=second, timestamp=2.0)
    store.add_channel_reading(channel=1, values_ma=third, timestamp=3.0)

    history = store.get_history(channel=1, led_index=1)
    assert [point.timestamp for point in history] == [2.0, 3.0]
    assert [point.value_ma for point in history] == [11.0, 21.0]


def test_set_max_points_keeps_recent_history() -> None:
    store = DataStore(max_points=3)
    first = tuple(float(index) for index in range(1, 17))
    second = tuple(float(index + 10) for index in range(1, 17))
    third = tuple(float(index + 20) for index in range(1, 17))

    store.add_channel_reading(channel=1, values_ma=first, timestamp=1.0)
    store.add_channel_reading(channel=1, values_ma=second, timestamp=2.0)
    store.add_channel_reading(channel=1, values_ma=third, timestamp=3.0)
    store.set_max_points(2)

    history = store.get_history(channel=1, led_index=1)
    assert [point.timestamp for point in history] == [2.0, 3.0]
    assert [point.value_ma for point in history] == [11.0, 21.0]


def test_clear_history_only_clears_selected_led() -> None:
    store = DataStore()
    values = tuple(float(index) for index in range(1, 17))

    store.add_channel_reading(channel=1, values_ma=values, timestamp=1.0)
    store.clear_history(channel=1, led_index=1)

    assert store.get_history(channel=1, led_index=1) == []
    assert len(store.get_history(channel=1, led_index=2)) == 1
