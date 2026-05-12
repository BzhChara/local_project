from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from time import time


DEFAULT_MAX_POINTS = 600


@dataclass(frozen=True)
class HistoryPoint:
    """单个历史采样点，timestamp 使用系统 Unix 时间戳。"""

    timestamp: float
    value_ma: float


class DataStore:
    """保存通道最新值和每个 LED 的内存历史曲线数据。

    当前阶段只做运行时内存缓存，不落盘保存。为避免长时间运行无限占用内存，
    每个通道/LED 默认最多保留最近 600 个点。
    """

    def __init__(self, max_points: int = DEFAULT_MAX_POINTS) -> None:
        self.max_points = max_points
        self.latest_readings: dict[int, tuple[float, ...]] = {}
        self._history: defaultdict[tuple[int, int], deque[HistoryPoint]] = defaultdict(
            lambda: deque(maxlen=self.max_points)
        )

    def add_channel_reading(
        self,
        channel: int,
        values_ma: tuple[float, ...],
        timestamp: float | None = None,
    ) -> None:
        """写入一个通道的 16 路数据，并同步追加每路 LED 的历史点。"""
        sample_time = time() if timestamp is None else timestamp
        self.latest_readings[channel] = values_ma
        for led_index, value_ma in enumerate(values_ma, start=1):
            self._history[(channel, led_index)].append(HistoryPoint(sample_time, value_ma))

    def get_latest(self, channel: int) -> tuple[float, ...] | None:
        """获取指定通道的最新 16 路数据。"""
        return self.latest_readings.get(channel)

    def get_history(self, channel: int, led_index: int) -> list[HistoryPoint]:
        """获取指定通道/LED 的历史点列表，供曲线刷新使用。"""
        return list(self._history[(channel, led_index)])
