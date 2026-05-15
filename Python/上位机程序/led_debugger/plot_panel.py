from __future__ import annotations

from time import time

import pyqtgraph as pg
from PySide6.QtWidgets import QLabel, QVBoxLayout

from led_debugger.data_store import HistoryPoint
from led_debugger.settings import AppSettings
from led_debugger.widgets import AcrylicPanel


AUTO_Y_RANGE_PADDING_RATIO = 0.15
AUTO_Y_SINGLE_VALUE_RATIO = 0.0025
AUTO_Y_MIN_PADDING = 0.00002
CURVE_STEP_MODE = "left"
CURVE_GRID_ALPHA = 0.12
CURVE_AXIS_COLOR = "#94a3b8"
CURVE_TICK_COLOR = "#475569"
CURVE_LINE_COLOR = "#1d4ed8"


class CurvePlotPanel(AcrylicPanel):
    """右侧 LED 电流历史曲线面板。"""

    def __init__(self, settings: AppSettings) -> None:
        super().__init__(radius=20, top_alpha=118, bottom_alpha=88)
        self.settings = settings
        self.setObjectName("curvePanel")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        self.curve_title = QLabel("LED电流 1 实时曲线")
        self.curve_title.setObjectName("curveTitle")

        self.plot_widget = pg.PlotWidget(axisItems={"bottom": pg.DateAxisItem(orientation="bottom")})
        self.plot_widget.setObjectName("curvePlot")
        self.plot_widget.setMinimumHeight(260)
        self.plot_widget.setBackground(None)
        self.plot_widget.setStyleSheet("background: transparent; border: none;")
        self.plot_widget.viewport().setStyleSheet("background: transparent;")
        self.plot_widget.showGrid(x=self.settings.curve_show_grid, y=self.settings.curve_show_grid, alpha=CURVE_GRID_ALPHA)
        self.plot_widget.setLabel("left", "电流值", units="mA", color=CURVE_TICK_COLOR)
        self.plot_widget.setLabel("bottom", "时间", color=CURVE_TICK_COLOR)
        self.plot_widget.setMouseEnabled(x=False, y=False)
        self.plot_widget.getPlotItem().setMenuEnabled(False)
        self.plot_widget.getPlotItem().hideButtons()
        self.plot_widget.getAxis("left").setPen(pg.mkPen(CURVE_AXIS_COLOR))
        self.plot_widget.getAxis("bottom").setPen(pg.mkPen(CURVE_AXIS_COLOR))
        self.plot_widget.getAxis("left").setTextPen(pg.mkPen(CURVE_TICK_COLOR))
        self.plot_widget.getAxis("bottom").setTextPen(pg.mkPen(CURVE_TICK_COLOR))
        self.plot_widget.getAxis("left").setStyle(tickTextOffset=8)
        self.plot_widget.getAxis("bottom").setStyle(tickTextOffset=8)
        self.plot_curve = self.plot_widget.plot(
            [],
            [],
            pen=self._curve_pen(),
            antialias=False,
            stepMode=CURVE_STEP_MODE,
        )
        self._reset_curve_axes()

        layout.addWidget(self.curve_title)
        layout.addWidget(self.plot_widget, 1)

    def set_selected_led(self, led_index: int) -> None:
        """更新曲线标题中的 LED 序号。"""
        self.curve_title.setText(f"LED电流 {led_index} 实时曲线")

    def apply_settings(self, settings: AppSettings) -> None:
        """应用曲线显示设置。"""
        self.settings = settings
        self.plot_widget.showGrid(
            x=self.settings.curve_show_grid,
            y=self.settings.curve_show_grid,
            alpha=CURVE_GRID_ALPHA,
        )
        self.plot_curve.setPen(self._curve_pen())
        x_values, y_values = self.plot_curve.getData()
        if x_values is None or y_values is None:
            self.plot_curve.setData([], [], antialias=False, stepMode=CURVE_STEP_MODE)
            self._apply_curve_y_range([])
        else:
            self.plot_curve.setData(x_values, y_values, antialias=False, stepMode=CURVE_STEP_MODE)
            self._apply_curve_y_range(list(y_values))

    def set_history(self, points: list[HistoryPoint]) -> None:
        """刷新当前通道、当前 LED 的实时历史曲线。"""
        if not points:
            self.plot_curve.setData([], [], antialias=False, stepMode=CURVE_STEP_MODE)
            self._reset_curve_axes()
            return

        latest_timestamp = points[-1].timestamp
        window_start = latest_timestamp - self.settings.curve_visible_seconds
        visible_points = [point for point in points if point.timestamp >= window_start]
        x_values = [point.timestamp for point in visible_points]
        y_values = [point.value_ma for point in visible_points]
        self.plot_curve.setData(x_values, y_values, antialias=False, stepMode=CURVE_STEP_MODE)
        self.plot_widget.setXRange(window_start, latest_timestamp, padding=0)
        self._apply_curve_y_range(y_values)

    def _curve_pen(self):
        """返回整数像素线宽的曲线画笔，减少亚像素抗锯齿伪影。"""
        return pg.mkPen(CURVE_LINE_COLOR, width=max(1, round(self.settings.curve_line_width)))

    def _reset_curve_axes(self) -> None:
        """无数据时显示最近一段真实时间，避免坐标轴落到 Unix 0 秒。"""
        latest_timestamp = time()
        window_start = latest_timestamp - self.settings.curve_visible_seconds
        self.plot_widget.setXRange(window_start, latest_timestamp, padding=0)
        self._apply_curve_y_range([])

    def _apply_curve_y_range(self, y_values: list[float] | None = None) -> None:
        """根据设置应用曲线 Y 轴范围。"""
        if not self.settings.y_axis_auto_scale:
            self.plot_widget.getPlotItem().enableAutoRange(axis="y", enable=False)
            self.plot_widget.setYRange(self.settings.y_axis_min, self.settings.y_axis_max, padding=0)
            return

        self.plot_widget.getPlotItem().enableAutoRange(axis="y", enable=False)

        if not y_values:
            self.plot_widget.setYRange(0.0, 1.0, padding=0)
            return

        minimum = min(y_values)
        maximum = max(y_values)
        if minimum == maximum:
            padding = self._auto_y_single_value_padding(minimum)
        else:
            padding = max((maximum - minimum) * AUTO_Y_RANGE_PADDING_RATIO, AUTO_Y_MIN_PADDING)
        self.plot_widget.setYRange(minimum - padding, maximum + padding, padding=0)

    def _auto_y_single_value_padding(self, value: float) -> float:
        """稳定电流值时给很小边距，避免 20mA 被放大到十几到二十几。"""
        return max(abs(value) * AUTO_Y_SINGLE_VALUE_RATIO, AUTO_Y_MIN_PADDING)
