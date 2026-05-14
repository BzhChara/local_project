from __future__ import annotations

import sys

import pyqtgraph as pg
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import (
    QColor,
    QFont,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QRadialGradient,
)
from PySide6.QtWidgets import (
    QApplication,
    QFrame,
    QGraphicsDropShadowEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
)

from led_debugger.data_store import DataStore
from led_debugger.serial_worker import SerialPoller, list_serial_ports


def apply_shadow(widget: QWidget, blur_radius: float, y_offset: float, alpha: int) -> QGraphicsDropShadowEffect:
    """给控件添加轻量阴影，用来模拟 Fluent/Acrylic 的悬浮层次。"""
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setBlurRadius(blur_radius)
    shadow.setOffset(0, y_offset)
    shadow.setColor(QColor(15, 23, 42, alpha))
    widget.setGraphicsEffect(shadow)
    return shadow


class AcrylicWindowBackground(QWidget):
    """带浅色渐变和雾化光斑的窗口背景。

    Qt 不能直接对桌面背景做系统级 Acrylic Blur，这里通过可透出的底层光斑、
    半透明玻璃层和高光边框来模拟 Windows 11 白色亚克力质感。
    """

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = self.rect()
        base = QLinearGradient(rect.topLeft(), rect.bottomRight())
        base.setColorAt(0.0, QColor("#f8fafc"))
        base.setColorAt(1.0, QColor("#e8eef7"))
        painter.fillRect(rect, base)

        # 大面积低透明度光斑提供“玻璃层背后有内容”的视觉基础。
        self._paint_blob(painter, rect.width() * 0.10, rect.height() * 0.12, rect.width() * 0.34, QColor(96, 165, 250, 46))
        self._paint_blob(painter, rect.width() * 0.82, rect.height() * 0.14, rect.width() * 0.30, QColor(168, 85, 247, 26))
        self._paint_blob(painter, rect.width() * 0.18, rect.height() * 0.86, rect.width() * 0.34, QColor(45, 212, 191, 31))
        self._paint_blob(painter, rect.width() * 0.56, rect.height() * 0.52, rect.width() * 0.42, QColor(96, 165, 250, 24))

        # 极轻的噪声点避免背景过于“塑料白”，但透明度很低，不影响阅读。
        painter.setPen(QPen(QColor(255, 255, 255, 22), 1))
        step = 18
        for x in range(0, rect.width(), step):
            offset = (x // step) % 3 * 6
            for y in range(offset, rect.height(), step * 2):
                painter.drawPoint(x, y)

        super().paintEvent(event)

    def _paint_blob(self, painter: QPainter, cx: float, cy: float, radius: float, color: QColor) -> None:
        """绘制一个从中心向外透明的柔和光斑。"""
        gradient = QRadialGradient(cx, cy, radius)
        gradient.setColorAt(0.0, color)
        gradient.setColorAt(0.55, QColor(color.red(), color.green(), color.blue(), int(color.alpha() * 0.42)))
        gradient.setColorAt(1.0, QColor(color.red(), color.green(), color.blue(), 0))
        painter.setPen(Qt.NoPen)
        painter.setBrush(gradient)
        painter.drawEllipse(QRectF(cx - radius, cy - radius, radius * 2, radius * 2))


class AcrylicPanel(QFrame):
    """可复用的伪亚克力面板。

    通过半透明渐变填充、顶部高光和白色细边框模拟玻璃厚度感。
    """

    def __init__(self, radius: int = 24, top_alpha: int = 138, bottom_alpha: int = 108, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.radius = radius
        self.top_alpha = top_alpha
        self.bottom_alpha = bottom_alpha
        self.setAttribute(Qt.WA_TranslucentBackground)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        path = QPainterPath()
        path.addRoundedRect(rect, self.radius, self.radius)

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, self.top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, self.bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        highlight = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        highlight.setColorAt(0.0, QColor(255, 255, 255, 92))
        highlight.setColorAt(0.42, QColor(255, 255, 255, 18))
        highlight.setColorAt(1.0, QColor(255, 255, 255, 0))
        painter.setBrush(highlight)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(QColor(255, 255, 255, 190), 1.0))
        painter.drawPath(path)

        super().paintEvent(event)


class AcrylicButton(QPushButton):
    """玻璃质感箭头按钮。"""

    def __init__(self, text: str) -> None:
        super().__init__(text)
        self._hovered = False
        self.setCursor(Qt.PointingHandCursor)
        self.setAttribute(Qt.WA_TranslucentBackground)

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self.update()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self.update()
        super().leaveEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        path = QPainterPath()
        path.addRoundedRect(rect, 16, 16)

        pressed = self.isDown()
        top_alpha = 168 if self._hovered else 132
        bottom_alpha = 126 if self._hovered else 98
        if pressed:
            top_alpha -= 26
            bottom_alpha -= 20

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(QColor(255, 255, 255, 210 if self._hovered else 178), 1.0))
        painter.drawPath(path)

        font = QFont(self.font())
        font.setPointSize(24)
        font.setWeight(QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.setPen(QColor("#0f172a"))
        painter.drawText(rect, Qt.AlignCenter, self.text())


class AcrylicTextButton(QPushButton):
    """玻璃质感文本按钮。

    COM 选择框和连接按钮共用这一套绘制逻辑，避免不同系统控件
    混用后出现高度、圆角和边框不一致的问题。
    """

    def __init__(self, text: str) -> None:
        super().__init__(text)
        self._hovered = False
        self.setCursor(Qt.PointingHandCursor)
        self.setAttribute(Qt.WA_TranslucentBackground)

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self.update()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self.update()
        super().leaveEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        path = QPainterPath()
        path.addRoundedRect(rect, 12, 12)

        pressed = self.isDown()
        top_alpha = 176 if self._hovered else 154
        bottom_alpha = 132 if self._hovered else 112
        if pressed:
            top_alpha -= 22
            bottom_alpha -= 18

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        border_alpha = 230 if self._hovered else 194
        painter.setPen(QPen(QColor(203, 213, 225, border_alpha), 1.0))
        painter.drawPath(path)

        painter.setFont(self.font())
        painter.setPen(QColor("#0f172a"))
        painter.drawText(rect.adjusted(10, 0, -10, 0), Qt.AlignCenter, self.text())


class AcrylicSelectButton(AcrylicTextButton):
    """自定义 COM 选择按钮。

    不再使用 QComboBox，避免组合框主体和弹出层分别受系统默认样式影响。
    后续串口扫描时可通过 set_items 填充真实 COM 口列表。
    """

    def __init__(self, placeholder: str) -> None:
        super().__init__(placeholder)
        self._placeholder = placeholder
        self._items: list[str] = []
        self._popup: AcrylicPanel | None = None
        self._before_popup = None
        self.clicked.connect(self._show_popup)

    def currentText(self) -> str:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        return self.text()

    def set_items(self, items: list[str]) -> None:
        """更新可选 COM 口列表，列表为空时保留占位文本。"""
        self._items = items
        if self.text() not in items:
            self.setText(self._placeholder)

    def set_before_popup_callback(self, callback) -> None:
        """设置弹出前回调，用于实时刷新 COM 口列表。"""
        self._before_popup = callback

    def _show_popup(self) -> None:
        """显示主窗口内部浮层，避免系统弹出窗口自带深色阴影。"""
        if self._popup is not None:
            self._popup.deleteLater()
            self._popup = None
            return
        if self._before_popup is not None:
            self._before_popup()

        root = self.window().centralWidget()
        popup = AcrylicPanel(radius=12, top_alpha=176, bottom_alpha=138, parent=root)
        popup.setObjectName("selectPopup")
        apply_shadow(popup, blur_radius=18, y_offset=6, alpha=20)

        layout = QVBoxLayout(popup)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        if self._items:
            for item in self._items:
                option = QPushButton(item)
                option.setObjectName("selectOption")
                option.setCursor(Qt.PointingHandCursor)
                option.clicked.connect(lambda _checked=False, value=item: self._select_item(value))
                layout.addWidget(option)
        else:
            empty_label = QLabel("暂无可用COM口")
            empty_label.setObjectName("selectEmpty")
            empty_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(empty_label)

        popup.setFixedWidth(self.width())
        popup.adjustSize()
        popup.move(root.mapFromGlobal(self.mapToGlobal(self.rect().bottomLeft())))
        popup.raise_()
        self._popup = popup
        popup.show()

    def _select_item(self, value: str) -> None:
        """选中一个 COM 口并关闭弹出层。"""
        self.setText(value)
        if self._popup is not None:
            self._popup.deleteLater()
            self._popup = None


class LedCard(QPushButton):
    """单个 LED 电流显示卡片。

    使用按钮控件便于后续点击后切换曲线。当前步骤先展示静态占位值，
    后续串口数据接入后只需要调用 set_current_value 更新显示。
    """

    def __init__(self, led_index: int) -> None:
        super().__init__()
        self.led_index = led_index
        self.setObjectName("ledCard")
        self.setCursor(Qt.PointingHandCursor)
        self.setCheckable(True)
        self.setMinimumSize(132, 104)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._hovered = False
        self._value_text = "--.-- mA"
        self._shadow = apply_shadow(self, blur_radius=18, y_offset=4, alpha=22)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.toggled.connect(lambda _checked: self._refresh_shadow())
        self.set_current_value(None)
        self._refresh_shadow()

    def set_current_value(self, value_ma: float | None) -> None:
        """更新卡片内显示的电流值。"""
        self._value_text = "--.-- mA" if value_ma is None else f"{value_ma:.3f} mA"
        self.setText(f"LED电流 {self.led_index}\n{self._value_text}")
        self.update()

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self._refresh_shadow()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self._refresh_shadow()
        super().leaveEvent(event)

    def _refresh_shadow(self) -> None:
        """根据普通、悬停、选中状态调整卡片阴影。"""
        if self.isChecked():
            blur_radius, y_offset, alpha = 30, 9, 46
        elif self._hovered:
            blur_radius, y_offset, alpha = 26, 7, 34
        else:
            blur_radius, y_offset, alpha = 22, 5, 24
        self._shadow.setBlurRadius(blur_radius)
        self._shadow.setOffset(0, y_offset)
        if self.isChecked():
            self._shadow.setColor(QColor(59, 130, 246, alpha))
        else:
            self._shadow.setColor(QColor(15, 23, 42, alpha))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        path = QPainterPath()
        path.addRoundedRect(rect, 18, 18)

        checked = self.isChecked()
        pressed = self.isDown()
        if checked:
            top_color = QColor(219, 234, 254, 150)
            bottom_color = QColor(219, 234, 254, 112)
            border_color = QColor(59, 130, 246, 218)
        else:
            top_alpha = 154 if self._hovered else 128
            bottom_alpha = 124 if self._hovered else 92
            if pressed:
                top_alpha -= 18
                bottom_alpha -= 16
            top_color = QColor(255, 255, 255, top_alpha)
            bottom_color = QColor(255, 255, 255, bottom_alpha)
            border_color = QColor(255, 255, 255, 230 if self._hovered else 178)

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, top_color)
        fill.setColorAt(1.0, bottom_color)
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        # 顶部高光让卡片看起来有一层玻璃厚度。
        highlight = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        highlight.setColorAt(0.0, QColor(255, 255, 255, 108))
        highlight.setColorAt(0.38, QColor(255, 255, 255, 18))
        highlight.setColorAt(1.0, QColor(255, 255, 255, 0))
        painter.setBrush(highlight)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(border_color, 1.5 if checked else 1.0))
        painter.drawPath(path)

        title_font = QFont(self.font())
        title_font.setPointSize(11)
        title_font.setWeight(QFont.Weight.Normal)
        value_font = QFont(self.font())
        value_font.setPointSize(15)
        value_font.setWeight(QFont.Weight.Normal)

        painter.setPen(QColor("#0f172a"))
        painter.setFont(title_font)
        painter.drawText(rect.adjusted(12, 16, -12, -rect.height() * 0.48), Qt.AlignHCenter | Qt.AlignTop, f"LED电流 {self.led_index}")

        painter.setFont(value_font)
        painter.drawText(rect.adjusted(12, rect.height() * 0.33, -12, -12), Qt.AlignCenter, self._value_text)


class MainWindow(QMainWindow):
    """主窗口。

    当前步骤实现完整静态布局。串口轮询和实时数据更新会在后续步骤接入，
    这样界面结构可以先被单独 review。
    """

    def __init__(self) -> None:
        super().__init__()
        self.current_channel = 1
        self.selected_led = 1
        self.led_cards: list[LedCard] = []
        self.available_ports: list[str] = []
        self.data_store = DataStore()
        self.serial_poller: SerialPoller | None = None
        self.setWindowTitle("LED亮度检测上位机")
        self.resize(1080, 720)
        self._build_ui()
        self._select_led(1)
        self._refresh_ports()

    def _build_ui(self) -> None:
        # 根容器自绘渐变和光斑，为半透明玻璃面板提供可透出的背景内容。
        root = AcrylicWindowBackground()
        root.setObjectName("root")
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(24, 20, 24, 20)
        root_layout.setSpacing(18)

        # 顶部工具栏：后续会填充真实 COM 口列表，并接入连接/断开逻辑。
        top_bar = AcrylicPanel(radius=22, top_alpha=118, bottom_alpha=86)
        top_bar.setObjectName("topBar")
        toolbar = QHBoxLayout(top_bar)
        toolbar.setContentsMargins(14, 10, 14, 10)
        toolbar.setSpacing(10)

        self.port_combo = AcrylicSelectButton("选择COM口")
        self.port_combo.setFixedSize(150, 38)
        self.port_combo.set_before_popup_callback(self._refresh_ports)

        self.connect_button = AcrylicTextButton("连接")
        self.connect_button.setObjectName("connectButton")
        self.connect_button.setFixedSize(92, 38)
        self.connect_button.clicked.connect(self._toggle_connection)

        self.status_label = QLabel("状态：未连接")
        self.status_label.setObjectName("statusLabel")

        toolbar.addWidget(self.port_combo)
        toolbar.addWidget(self.connect_button)
        toolbar.addWidget(self.status_label)
        toolbar.addItem(QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))

        apply_shadow(top_bar, blur_radius=30, y_offset=10, alpha=24)
        root_layout.addWidget(top_bar)

        # 主内容面板：左侧为 4x4 LED 电流卡片，右侧为选中 LED 的曲线区域。
        content = AcrylicPanel(radius=26, top_alpha=142, bottom_alpha=104)
        content.setObjectName("contentPanel")
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(20, 20, 20, 20)
        content_layout.setSpacing(18)

        work_area = QHBoxLayout()
        work_area.setSpacing(18)

        cards_panel = QFrame()
        cards_panel.setObjectName("cardsPanel")
        cards_panel_layout = QVBoxLayout(cards_panel)
        cards_panel_layout.setContentsMargins(0, 0, 0, 0)
        cards_panel_layout.setSpacing(14)

        self.channel_label = QLabel("通道 1")
        self.channel_label.setObjectName("channelTitle")
        self.channel_label.setAlignment(Qt.AlignCenter)

        cards_grid = QWidget()
        cards_grid.setObjectName("cardsGrid")
        cards_layout = QGridLayout(cards_grid)
        cards_layout.setContentsMargins(0, 0, 0, 0)
        cards_layout.setSpacing(12)

        for index in range(1, 17):
            card = LedCard(index)
            card.clicked.connect(lambda _checked=False, led=index: self._select_led(led))
            self.led_cards.append(card)
            row = (index - 1) // 4
            column = (index - 1) % 4
            cards_layout.addWidget(card, row, column)

        cards_panel_layout.addWidget(self.channel_label)
        cards_panel_layout.addWidget(cards_grid, 1)

        curve_panel = AcrylicPanel(radius=20, top_alpha=118, bottom_alpha=88)
        curve_panel.setObjectName("curvePanel")
        curve_layout = QVBoxLayout(curve_panel)
        curve_layout.setContentsMargins(18, 16, 18, 16)
        curve_layout.setSpacing(12)

        self.curve_title = QLabel("LED电流 1 实时曲线")
        self.curve_title.setObjectName("curveTitle")

        self.plot_widget = pg.PlotWidget(axisItems={"bottom": pg.DateAxisItem(orientation="bottom")})
        self.plot_widget.setObjectName("curvePlot")
        self.plot_widget.setMinimumHeight(260)
        self.plot_widget.setBackground((255, 255, 255, 24))
        self.plot_widget.showGrid(x=True, y=True, alpha=0.22)
        self.plot_widget.setLabel("left", "电流值", units="mA")
        self.plot_widget.setLabel("bottom", "时间")
        self.plot_widget.getAxis("left").setPen(pg.mkPen("#94a3b8"))
        self.plot_widget.getAxis("bottom").setPen(pg.mkPen("#94a3b8"))
        self.plot_widget.getAxis("left").setTextPen(pg.mkPen("#64748b"))
        self.plot_widget.getAxis("bottom").setTextPen(pg.mkPen("#64748b"))
        self.plot_curve = self.plot_widget.plot([], [], pen=pg.mkPen("#3b82f6", width=2))

        curve_layout.addWidget(self.curve_title)
        curve_layout.addWidget(self.plot_widget, 1)

        apply_shadow(content, blur_radius=44, y_offset=16, alpha=28)
        apply_shadow(curve_panel, blur_radius=30, y_offset=10, alpha=24)
        work_area.addWidget(cards_panel, 3)
        work_area.addWidget(curve_panel, 2)
        content_layout.addLayout(work_area, 1)

        # 通道切换按钮固定在面板底部左右两侧，符合用户提出的切换方式。
        nav_layout = QHBoxLayout()
        self.prev_button = AcrylicButton("‹")
        self.prev_button.setObjectName("navButton")
        self.prev_button.setFixedSize(56, 44)
        self.next_button = AcrylicButton("›")
        self.next_button.setObjectName("navButton")
        self.next_button.setFixedSize(56, 44)
        self.prev_button.clicked.connect(lambda: self._change_channel(-1))
        self.next_button.clicked.connect(lambda: self._change_channel(1))
        apply_shadow(self.prev_button, blur_radius=22, y_offset=7, alpha=28)
        apply_shadow(self.next_button, blur_radius=22, y_offset=7, alpha=28)
        nav_layout.addWidget(self.prev_button, 0, Qt.AlignLeft)
        nav_layout.addItem(QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
        nav_layout.addWidget(self.next_button, 0, Qt.AlignRight)
        content_layout.addLayout(nav_layout)

        root_layout.addWidget(content, 1)
        self.setCentralWidget(root)
        # 样式集中在当前窗口内，后续如果界面复杂再拆分为独立 qss 文件。
        self.setStyleSheet(
            """
            #root {
                background: transparent;
                color: #0f172a;
                font-family: "Microsoft YaHei UI", "Segoe UI", sans-serif;
                font-size: 14px;
            }
            #connectButton {
                background: transparent;
                border: none;
                color: #0f172a;
                padding: 0;
            }
            #selectOption {
                min-height: 30px;
                border: none;
                border-radius: 8px;
                background: transparent;
                color: #0f172a;
                text-align: center;
            }
            #selectOption:hover {
                background: rgba(219, 234, 254, 0.86);
            }
            #selectEmpty {
                min-height: 30px;
                color: #475569;
                background: transparent;
            }
            #statusLabel {
                color: #64748b;
            }
            #channelTitle {
                color: #0f172a;
                font-size: 24px;
                font-weight: 600;
            }
            #cardsPanel {
                background: transparent;
            }
            #cardsGrid {
                background: transparent;
            }
            #ledCard {
                color: #0f172a;
                font-size: 15px;
                font-weight: 400;
                background: transparent;
                border: none;
            }
            #curveTitle {
                color: #0f172a;
                font-size: 17px;
                font-weight: 600;
            }
            #curvePlaceholder {
                color: #64748b;
                background: rgba(255, 255, 255, 0.22);
                border: 1px dashed rgba(148, 163, 184, 0.35);
                border-radius: 14px;
                font-size: 14px;
            }
            #navButton {
                background: transparent;
                border: none;
                padding: 0;
            }
            """
        )

    def _change_channel(self, offset: int) -> None:
        """循环切换 1-16 通道，并刷新界面标题。"""
        self.current_channel = ((self.current_channel - 1 + offset) % 16) + 1
        self.channel_label.setText(f"通道 {self.current_channel}")
        self._update_cards_for_current_channel()
        self._update_curve()

    def _select_led(self, led_index: int) -> None:
        """选中一个 LED 卡片，并更新曲线区域标题。"""
        self.selected_led = led_index
        for card in self.led_cards:
            card.setChecked(card.led_index == led_index)
        self.curve_title.setText(f"LED电流 {led_index} 实时曲线")
        self._update_curve()

    def _refresh_ports(self) -> None:
        """刷新可用 COM 口列表。"""
        self.available_ports = list_serial_ports()
        self.port_combo.set_items(self.available_ports)
        if self.serial_poller is None:
            if self.available_ports:
                self.status_label.setText(f"状态：未连接，发现 {len(self.available_ports)} 个COM口")
            else:
                self.status_label.setText("状态：未连接，未发现COM口")

    def _toggle_connection(self) -> None:
        """连接或断开串口轮询线程。"""
        if self.serial_poller is not None:
            self._disconnect_serial()
            return

        self._refresh_ports()
        port_name = self.port_combo.currentText()
        if port_name not in self.available_ports:
            self.status_label.setText("状态：请选择COM口")
            return

        self.serial_poller = SerialPoller(port_name, self)
        self.serial_poller.reading_received.connect(self._handle_reading_received)
        self.serial_poller.status_changed.connect(self.status_label.setText)
        self.serial_poller.error_occurred.connect(self._handle_serial_error)
        self.serial_poller.finished.connect(self._handle_serial_finished)
        self.connect_button.setText("断开")
        self.status_label.setText(f"状态：连接中 {port_name}")
        self.serial_poller.start()

    def _disconnect_serial(self) -> None:
        """停止串口线程，按钮状态在 finished 回调中恢复。"""
        if self.serial_poller is None:
            return
        self.status_label.setText("状态：正在断开...")
        self.connect_button.setEnabled(False)
        self.serial_poller.stop()

    def _handle_reading_received(self, channel: int, currents_ma: tuple[float, ...]) -> None:
        """接收后台线程推送的最新通道数据，并刷新当前页面。"""
        self.data_store.add_channel_reading(channel, currents_ma)
        if channel == self.current_channel:
            self._update_cards_for_current_channel()
            self._update_curve()

    def _handle_serial_error(self, message: str) -> None:
        """显示串口错误；轮询线程仍会继续尝试后续通道。"""
        self.status_label.setText(f"状态：{message}")

    def _handle_serial_finished(self) -> None:
        """串口线程退出后恢复连接按钮状态。"""
        self.serial_poller = None
        self.connect_button.setText("连接")
        self.connect_button.setEnabled(True)

    def _update_cards_for_current_channel(self) -> None:
        """用当前通道最新值刷新 16 个 LED 卡片。"""
        values = self.data_store.get_latest(self.current_channel)
        for index, card in enumerate(self.led_cards):
            card.set_current_value(values[index] if values is not None else None)

    def _update_curve(self) -> None:
        """刷新当前通道、当前 LED 的实时历史曲线。"""
        points = self.data_store.get_history(self.current_channel, self.selected_led)
        if not points:
            self.plot_curve.setData([], [])
            return
        x_values = [point.timestamp for point in points]
        y_values = [point.value_ma for point in points]
        self.plot_curve.setData(x_values, y_values)

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        """关闭窗口前停止后台串口线程，避免进程残留。"""
        if self.serial_poller is not None:
            self.serial_poller.stop()
            self.serial_poller.wait(1000)
        super().closeEvent(event)


def run() -> int:
    """创建 Qt 应用并显示主窗口。"""
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()
