from __future__ import annotations

import ctypes
import os
import sys
import time

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QIcon, QLinearGradient, QPainter, QPen, QRadialGradient, QTextCursor
from PySide6.QtWidgets import (
    QApplication,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QSizePolicy,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from app_paths import get_app_dir
from config_loader import DEFAULT_OPTIONS, get_config_path, load_config, parse_bool
from filter_logic import WAVELENGTH_COLORS, channel_indices
from spectro_ui.serial_worker import SerialTestWorker, TestConfig, list_serial_port_names
from spectro_ui.widgets import AcrylicPanel, AcrylicSelectButton, AcrylicTextButton, ChannelCard, ElidedLabel, apply_shadow


APP_ICON_FILE = os.path.join("assets", "app_icon_spectro_transparent.png")


def load_app_icon() -> QIcon:
    for base_dir in (get_app_dir(), getattr(sys, "_MEIPASS", None)):
        if not base_dir:
            continue
        icon_path = os.path.join(base_dir, APP_ICON_FILE)
        if os.path.exists(icon_path):
            return QIcon(icon_path)
    return QIcon()


def _colorref(red: int, green: int, blue: int) -> int:
    return red | (green << 8) | (blue << 16)


def apply_native_title_bar_style(window: QWidget) -> None:
    if sys.platform != "win32":
        return

    try:
        hwnd = int(window.winId())
        caption_color = ctypes.c_int(_colorref(248, 250, 252))
        text_color = ctypes.c_int(_colorref(15, 23, 42))
        border_color = ctypes.c_int(_colorref(226, 232, 240))
        dwm = ctypes.windll.dwmapi
        dwm.DwmSetWindowAttribute(hwnd, 35, ctypes.byref(caption_color), ctypes.sizeof(caption_color))
        dwm.DwmSetWindowAttribute(hwnd, 36, ctypes.byref(text_color), ctypes.sizeof(text_color))
        dwm.DwmSetWindowAttribute(hwnd, 34, ctypes.byref(border_color), ctypes.sizeof(border_color))
    except Exception:
        return


class AcrylicWindowBackground(QWidget):
    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = self.rect()
        base = QLinearGradient(rect.topLeft(), rect.bottomRight())
        base.setColorAt(0.0, QColor("#f8fafc"))
        base.setColorAt(1.0, QColor("#e8eef7"))
        painter.fillRect(rect, base)

        self._paint_blob(painter, rect.width() * 0.10, rect.height() * 0.12, rect.width() * 0.34, QColor(96, 165, 250, 46))
        self._paint_blob(painter, rect.width() * 0.82, rect.height() * 0.14, rect.width() * 0.30, QColor(168, 85, 247, 26))
        self._paint_blob(painter, rect.width() * 0.18, rect.height() * 0.86, rect.width() * 0.34, QColor(45, 212, 191, 31))
        self._paint_blob(painter, rect.width() * 0.56, rect.height() * 0.52, rect.width() * 0.42, QColor(96, 165, 250, 24))

        painter.setPen(QPen(QColor(255, 255, 255, 22), 1))
        step = 18
        for x in range(0, rect.width(), step):
            offset = (x // step) % 3 * 6
            for y in range(offset, rect.height(), step * 2):
                painter.drawPoint(x, y)

        super().paintEvent(event)

    def _paint_blob(self, painter: QPainter, cx: float, cy: float, radius: float, color: QColor) -> None:
        gradient = QRadialGradient(cx, cy, radius)
        gradient.setColorAt(0.0, color)
        gradient.setColorAt(0.55, QColor(color.red(), color.green(), color.blue(), int(color.alpha() * 0.42)))
        gradient.setColorAt(1.0, QColor(color.red(), color.green(), color.blue(), 0))
        painter.setPen(Qt.NoPen)
        painter.setBrush(gradient)
        painter.drawEllipse(QRectF(cx - radius, cy - radius, radius * 2, radius * 2))


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("分光光度滤光片记录")
        self.setWindowIcon(load_app_icon())
        self.resize(1280, 760)
        self.worker: SerialTestWorker | None = None
        self.channel_cards: list[ChannelCard] = []
        self.selected_channel = 0
        self.active_channel_indices = list(range(24))
        self.config_path = get_config_path()
        self.initial_options = self._load_initial_options()

        self._build_ui()
        self._refresh_ports()
        self._apply_initial_options()
        self._set_running(False)
        apply_native_title_bar_style(self)

    def showEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        super().showEvent(event)
        apply_native_title_bar_style(self)

    def _build_ui(self) -> None:
        root = AcrylicWindowBackground()
        root.setObjectName("root")
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(24, 20, 24, 20)
        root_layout.setSpacing(18)

        top_bar = AcrylicPanel(radius=22, top_alpha=118, bottom_alpha=86)
        top_bar.setObjectName("topBar")
        toolbar = QHBoxLayout(top_bar)
        toolbar.setContentsMargins(14, 10, 14, 10)
        toolbar.setSpacing(10)

        self.port_combo = self._create_combo("选择COM口")
        self.port_combo.setFixedSize(150, 38)
        self.refresh_button = AcrylicTextButton("刷新串口")
        self.refresh_button.setFixedSize(92, 38)
        self.refresh_button.clicked.connect(self._refresh_ports)

        self.wavelength_combo = self._create_combo("选择波长")
        self.wavelength_combo.setFixedSize(150, 38)
        for wavelength in (410, 460, 520, 550, 590, 630):
            self.wavelength_combo.addItem(f"{wavelength}nm {WAVELENGTH_COLORS[wavelength]}光", wavelength)

        self.ratio_combo = self._create_combo("选择滤光片")
        self.ratio_combo.setFixedSize(132, 38)
        for ratio, label in ((0, "空气记录"), (10, "10% 滤光片"), (20, "20% 滤光片"), (30, "30% 滤光片")):
            self.ratio_combo.addItem(label, ratio)

        self.group_combo = self._create_combo("选择通道")
        self.group_combo.setFixedSize(132, 38)
        self.group_combo.addItem("CH1-CH24", 0)
        self.group_combo.addItem("CH1-CH12", 1)
        self.group_combo.addItem("CH13-CH24", 2)
        self.group_combo.set_selection_changed_callback(self._handle_group_changed)

        self.start_button = AcrylicTextButton("开始测试")
        self.start_button.setFixedSize(92, 38)
        self.start_button.clicked.connect(self._start_test)
        self.stop_button = AcrylicTextButton("停止")
        self.stop_button.setFixedSize(92, 38)
        self.stop_button.clicked.connect(self._stop_test)

        self.status_label = ElidedLabel("状态：未运行")
        self.status_label.setObjectName("statusLabel")
        self.status_label.setMinimumWidth(150)

        self.summary_label = QLabel("已写入 0 行，稳定 0/24")
        self.summary_label.setObjectName("summaryLabel")
        self.summary_label.setMinimumWidth(180)

        toolbar.addWidget(self.port_combo)
        toolbar.addWidget(self.refresh_button)
        toolbar.addWidget(self.wavelength_combo)
        toolbar.addWidget(self.ratio_combo)
        toolbar.addWidget(self.group_combo)
        toolbar.addWidget(self.start_button)
        toolbar.addWidget(self.stop_button)
        toolbar.addWidget(self.status_label, 1)
        toolbar.addWidget(self.summary_label)

        apply_shadow(top_bar, blur_radius=30, y_offset=10, alpha=24)
        root_layout.addWidget(top_bar)

        content = AcrylicPanel(radius=26, top_alpha=142, bottom_alpha=104)
        content.setObjectName("contentPanel")
        content_layout = QHBoxLayout(content)
        content_layout.setContentsMargins(20, 20, 20, 20)
        content_layout.setSpacing(18)

        cards_panel = QFrame()
        cards_panel.setObjectName("cardsPanel")
        cards_layout = QGridLayout(cards_panel)
        cards_layout.setContentsMargins(0, 0, 0, 0)
        cards_layout.setSpacing(12)

        for index in range(1, 25):
            card = ChannelCard(index)
            self.channel_cards.append(card)
            row = (index - 1) // 6
            column = (index - 1) % 6
            cards_layout.addWidget(card, row, column)

        log_panel = AcrylicPanel(radius=20, top_alpha=118, bottom_alpha=88)
        log_panel.setObjectName("logPanel")
        apply_shadow(log_panel, blur_radius=30, y_offset=10, alpha=24)
        log_layout = QVBoxLayout(log_panel)
        log_layout.setContentsMargins(18, 16, 18, 16)
        log_layout.setSpacing(12)

        log_title_layout = QHBoxLayout()
        log_title = QLabel("运行日志")
        log_title.setObjectName("sectionTitle")
        self.clear_log_button = AcrylicTextButton("清空")
        self.clear_log_button.clicked.connect(self._clear_log)
        log_title_layout.addWidget(log_title)
        log_title_layout.addStretch(1)
        log_title_layout.addWidget(self.clear_log_button)

        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setObjectName("logEdit")
        log_layout.addLayout(log_title_layout)
        log_layout.addWidget(self.log_edit, 1)

        content_layout.addWidget(cards_panel, 5)
        content_layout.addWidget(log_panel, 2)
        root_layout.addWidget(content, 1)

        self.setCentralWidget(root)
        self.setStyleSheet(
            """
            QWidget {
                background: transparent;
                color: #0f172a;
                font-family: "Microsoft YaHei UI", "Segoe UI", sans-serif;
                font-size: 13px;
            }
            #root {
                background: transparent;
            }
            #contentPanel {
                background: transparent;
            }
            #cardsPanel {
                background: transparent;
            }
            QLabel#statusLabel {
                color: #64748b;
            }
            QLabel#summaryLabel {
                color: #0f172a;
                font-weight: 600;
            }
            QLabel#sectionTitle {
                font-size: 17px;
                font-weight: 600;
            }
            QTextEdit#logEdit {
                border: 1px solid rgba(148, 163, 184, 0.36);
                border-radius: 12px;
                padding: 10px;
                background: rgba(255, 255, 255, 0.44);
                color: #0f172a;
                font-family: "Consolas", "Microsoft YaHei UI", monospace;
                font-size: 12px;
            }
            QTextEdit#logEdit QScrollBar:vertical {
                width: 8px;
                margin: 4px 0 4px 0;
                border: none;
                border-radius: 4px;
                background: rgba(226, 232, 240, 0.36);
            }
            QTextEdit#logEdit QScrollBar::handle:vertical {
                min-height: 44px;
                border-radius: 4px;
                background: rgba(100, 116, 139, 0.42);
            }
            QTextEdit#logEdit QScrollBar::handle:vertical:hover {
                background: rgba(71, 85, 105, 0.58);
            }
            QTextEdit#logEdit QScrollBar::handle:vertical:pressed {
                background: rgba(51, 65, 85, 0.65);
            }
            QTextEdit#logEdit QScrollBar::add-line:vertical,
            QTextEdit#logEdit QScrollBar::sub-line:vertical {
                height: 0;
                border: none;
                background: transparent;
            }
            QTextEdit#logEdit QScrollBar::add-page:vertical,
            QTextEdit#logEdit QScrollBar::sub-page:vertical {
                background: transparent;
            }
            """
        )

    def _create_combo(self, placeholder: str) -> AcrylicSelectButton:
        combo = AcrylicSelectButton(placeholder)
        combo.setFixedHeight(36)
        combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        return combo

    def _load_initial_options(self):
        values = dict(DEFAULT_OPTIONS)
        values["filter_ratio"] = 30
        config_values = load_config(self.config_path)
        for key, raw_value in config_values.items():
            if key in ("wavelength", "channel_group", "filter_ratio"):
                values[key] = self._safe_int(raw_value, values[key])
            elif key in ("ratio_tolerance", "air_tolerance"):
                values[key] = self._safe_float(raw_value, values[key])
            elif key in ("no_start", "keep_light"):
                values[key] = self._safe_bool(raw_value, values[key])
            elif key in values:
                values[key] = raw_value
        return values

    def _apply_initial_options(self) -> None:
        self._set_combo_by_data(self.wavelength_combo, self.initial_options["wavelength"])
        self._set_combo_by_data(self.ratio_combo, self.initial_options["filter_ratio"])
        self._set_combo_by_data(self.group_combo, self.initial_options["channel_group"])
        port = self.initial_options["port"]
        if port and self.port_combo.findText(port) < 0:
            self.port_combo.addItem(port, port)
        self._set_combo_by_text(self.port_combo, port)
        self._refresh_idle_summary()

    def _handle_group_changed(self, text: str) -> None:
        if self.worker is not None:
            return
        self._refresh_idle_summary()

    def _refresh_idle_summary(self) -> None:
        self.active_channel_indices = channel_indices(self.group_combo.currentData())
        self._set_summary(0, 0)

    def _refresh_ports(self) -> None:
        current = self.port_combo.currentText()
        ports = list_serial_port_names()
        self.port_combo.clear()
        for port in ports:
            self.port_combo.addItem(port, port)
        if current and self.port_combo.findText(current) < 0:
            self.port_combo.addItem(current, current)
        if ports or current:
            self._set_combo_by_text(self.port_combo, current or ports[0])
        self._append_log(f"已刷新串口: {', '.join(ports) if ports else '未发现串口'}")

    def _start_test(self) -> None:
        port = self.port_combo.currentText().strip()
        if not port:
            self._set_status("状态：请选择串口")
            return

        self.active_channel_indices = channel_indices(self.group_combo.currentData())
        self._reset_cards()
        self._set_summary(0, 0)
        config = TestConfig(
            port=port,
            wavelength=self.wavelength_combo.currentData(),
            channel_group=self.group_combo.currentData(),
            filter_ratio=self.ratio_combo.currentData(),
            output=self.initial_options["output"],
            stable_output=self.initial_options["stable_output"],
            ratio_tolerance=float(self.initial_options["ratio_tolerance"]),
            air_tolerance=float(self.initial_options["air_tolerance"]),
            no_start=bool(self.initial_options["no_start"]),
            keep_light=bool(self.initial_options["keep_light"]),
        )

        self.worker = SerialTestWorker(config, self)
        self.worker.status_changed.connect(self._set_status)
        self.worker.log_added.connect(self._append_log)
        self.worker.values_received.connect(self._handle_values)
        self.worker.baseline_progress_changed.connect(self._handle_baseline_progress)
        self.worker.baseline_value_changed.connect(self._handle_baseline_value)
        self.worker.channel_state_changed.connect(self._handle_channel_state)
        self.worker.stable_value_changed.connect(self._handle_stable_value)
        self.worker.row_status_changed.connect(self._handle_row_status)
        self.worker.output_paths_ready.connect(self._handle_output_paths)
        self.worker.test_finished.connect(self._handle_test_finished)

        self._set_running(True)
        self._append_log(
            f"开始测试: {port}, {config.wavelength}nm, 滤光片 {config.filter_ratio}%, "
            f"通道组 {self.group_combo.currentText()}"
        )
        self.worker.start()

    def _stop_test(self) -> None:
        if self.worker is None:
            return
        self._set_status("状态：正在停止...")
        self.worker.stop()
        self.stop_button.setEnabled(False)

    def _handle_test_finished(self) -> None:
        if self.worker is not None:
            self.worker.deleteLater()
            self.worker = None
        self._set_running(False)

    def _set_running(self, running: bool) -> None:
        self.start_button.setEnabled(not running)
        self.stop_button.setEnabled(running)
        self.port_combo.setEnabled(not running)
        self.refresh_button.setEnabled(not running)
        self.wavelength_combo.setEnabled(not running)
        self.ratio_combo.setEnabled(not running)
        self.group_combo.setEnabled(not running)

    def _handle_values(self, values: list[float]) -> None:
        for idx in self.active_channel_indices:
            value = values[idx]
            self.channel_cards[idx].set_value(value)

    def _handle_baseline_progress(self, done: int, total: int) -> None:
        self.summary_label.setText(f"空气基底 {done}/{total}")
        self.summary_label.setToolTip(f"空气基底 {done}/{total}，基底阶段不写入原始数据")
        self._set_status("状态：请保持空气状态，正在采集空气基底")

    def _handle_baseline_value(self, idx: int, value: float) -> None:
        self.channel_cards[idx].set_baseline(value)

    def _handle_channel_state(self, idx: int, state: str) -> None:
        self.channel_cards[idx].set_state(state)

    def _handle_stable_value(self, idx: int, value: float, ratio_percent: float) -> None:
        self.channel_cards[idx].set_stable(value)
        self._set_status(f"状态：CH{idx + 1} 稳定 {value:.6f}，比例 {ratio_percent:.2f}%")

    def _handle_row_status(self, rows: int, stable_count: int) -> None:
        self._set_summary(rows, stable_count)

    def _set_summary(self, rows: int, stable_count: int) -> None:
        total = len(self.active_channel_indices)
        summary = f"已写入 {rows} 行，稳定 {stable_count}/{total}"
        self.summary_label.setText(summary)
        self.summary_label.setToolTip(summary)

    def _handle_output_paths(self, raw_path: str, stable_path: str) -> None:
        self._append_log(f"原始数据文件: {raw_path}")
        if self.ratio_combo.currentData() != 0:
            self._append_log(f"稳定值文件: {stable_path}")

    def _reset_cards(self) -> None:
        for card in self.channel_cards:
            card.reset_measurement()
        self.selected_channel = 0

    def _set_status(self, message: str) -> None:
        short_message = self._short_status_message(message)
        self.status_label.setText(short_message)
        self.status_label.setToolTip(message)

    def _short_status_message(self, message: str) -> str:
        text = message.removeprefix("状态：")
        if "请选择串口" in text:
            return "状态：请选择串口"
        if "正在停止" in text:
            return "状态：正在停止"
        if "请保持空气状态" in text or "建立空气基底" in text:
            return "状态：基底采集中"
        if text.startswith("基底 "):
            return "状态：基底采集中"
        if "基底完成" in text:
            return "状态：请插入滤光片"
        if "空气记录模式" in text:
            return "状态：空气记录中"
        if "已连接" in text:
            return "状态：已连接"
        if "未运行" in text or "未连接" in text:
            return "状态：未运行"
        if "串口异常" in text:
            return "状态：串口异常"
        if "保存异常" in text:
            return "状态：保存异常"
        if "被占用" in text or "不能使用同一个路径" in text:
            return "状态：文件错误"
        if text.startswith("CH") and "稳定" in text:
            channel = text.split(" ", 1)[0]
            return f"状态：{channel} 稳定"
        return message if len(message) <= 12 else f"状态：{text[:8]}..."

    def _append_log(self, message: str) -> None:
        timestamp = time.strftime("%H:%M:%S")
        self.log_edit.append(f"[{timestamp}] {message}")
        self.log_edit.moveCursor(QTextCursor.MoveOperation.End)

    def _clear_log(self) -> None:
        self.log_edit.clear()

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.worker is not None:
            self.worker.stop()
            self.worker.wait(1200)
        super().closeEvent(event)

    def _set_combo_by_data(self, combo: AcrylicSelectButton, value) -> None:
        index = combo.findData(value)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _set_combo_by_text(self, combo: AcrylicSelectButton, text: str) -> None:
        index = combo.findText(text)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _safe_int(self, value, fallback):
        try:
            return int(value)
        except (TypeError, ValueError):
            return fallback

    def _safe_float(self, value, fallback):
        try:
            return float(value)
        except (TypeError, ValueError):
            return fallback

    def _safe_bool(self, value, fallback):
        try:
            return parse_bool(str(value))
        except ValueError:
            return fallback


def run() -> int:
    app = QApplication(sys.argv)
    app.setWindowIcon(load_app_icon())
    window = MainWindow()
    window.show()
    return app.exec()
