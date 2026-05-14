from __future__ import annotations

import sys
import ctypes
from dataclasses import replace
from pathlib import Path
from time import time

import pyqtgraph as pg
from PySide6.QtCore import QEasingCurve, QPropertyAnimation, QRectF, Qt, QTimer, Property
from PySide6.QtGui import (
    QColor,
    QDoubleValidator,
    QFont,
    QIcon,
    QIntValidator,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QRadialGradient,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QFrame,
    QGraphicsDropShadowEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from led_debugger.data_store import DataStore
from led_debugger.serial_worker import SerialPoller, list_serial_ports
from led_debugger.settings import (
    CURVE_VISIBLE_SECONDS_RANGE,
    CURVE_LINE_WIDTH_RANGE,
    DEFAULT_CURVE_VISIBLE_SECONDS,
    DEFAULT_CURVE_LINE_WIDTH,
    DEFAULT_CURVE_SHOW_GRID,
    DEFAULT_DISPLAY_DECIMALS,
    DEFAULT_MAX_POINTS,
    DEFAULT_READ_TIMEOUT_SECONDS,
    DEFAULT_ROUND_DELAY_SECONDS,
    DEFAULT_SHOW_SERIAL_ERRORS,
    DEFAULT_Y_AXIS_AUTO_SCALE,
    DEFAULT_Y_AXIS_MAX,
    DEFAULT_Y_AXIS_MIN,
    DISPLAY_DECIMALS_RANGE,
    MAX_POINTS_RANGE,
    READ_TIMEOUT_RANGE,
    ROUND_DELAY_RANGE,
    Y_AXIS_RANGE,
    AppSettings,
    load_settings,
    save_settings,
)


APP_ICON_PATH = Path(__file__).resolve().parent / "assets" / "app_icon.ico"
CURVE_VISIBLE_SECONDS = DEFAULT_CURVE_VISIBLE_SECONDS
AUTO_Y_RANGE_PADDING_RATIO = 0.15
AUTO_Y_SINGLE_VALUE_RATIO = 0.0025
AUTO_Y_MIN_PADDING = 0.00002
CURVE_STEP_MODE = "left"


def apply_shadow(widget: QWidget, blur_radius: float, y_offset: float, alpha: int) -> QGraphicsDropShadowEffect:
    """给控件添加轻量阴影，用来模拟 Fluent/Acrylic 的悬浮层次。"""
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setBlurRadius(blur_radius)
    shadow.setOffset(0, y_offset)
    shadow.setColor(QColor(15, 23, 42, alpha))
    widget.setGraphicsEffect(shadow)
    return shadow


def _colorref(red: int, green: int, blue: int) -> int:
    """Windows COLORREF 使用 0x00BBGGRR 字节顺序。"""
    return red | (green << 8) | (blue << 16)


def apply_native_title_bar_style(window: QMainWindow) -> None:
    """在支持的 Windows 版本上把原生标题栏设置为浅色。"""
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


def _lerp(start: float, end: float, progress: float) -> float:
    """线性插值，用于按钮动画的颜色和阴影过渡。"""
    return start + (end - start) * progress


class AnimatedButtonMixin:
    """给自绘按钮补充 Win 风格 hover/press 动画。"""

    def _init_button_animation(self) -> None:
        self._hover_progress = 0.0
        self._press_progress = 0.0

        self._hover_animation = QPropertyAnimation(self, b"hoverProgress", self)
        self._hover_animation.setDuration(130)
        self._hover_animation.setEasingCurve(QEasingCurve.Type.OutCubic)

        self._press_animation = QPropertyAnimation(self, b"pressProgress", self)
        self._press_animation.setDuration(105)
        self._press_animation.setEasingCurve(QEasingCurve.Type.OutCubic)

    def _get_hover_progress(self) -> float:
        return self._hover_progress

    def _set_hover_progress(self, value: float) -> None:
        self._hover_progress = value
        self._update_animation_frame()

    hoverProgress = Property(float, _get_hover_progress, _set_hover_progress)

    def _get_press_progress(self) -> float:
        return self._press_progress

    def _set_press_progress(self, value: float) -> None:
        self._press_progress = value
        self._update_animation_frame()

    pressProgress = Property(float, _get_press_progress, _set_press_progress)

    def _animate_hover(self, target: float) -> None:
        self._hover_animation.stop()
        self._hover_animation.setStartValue(self._hover_progress)
        self._hover_animation.setEndValue(target)
        self._hover_animation.start()

    def _animate_press(self, target: float, duration: int | None = None) -> None:
        self._press_animation.stop()
        if duration is not None:
            self._press_animation.setDuration(duration)
        self._press_animation.setStartValue(self._press_progress)
        self._press_animation.setEndValue(target)
        self._press_animation.start()

    def _animated_rect(self, rect: QRectF, scale_ratio: float = 0.006) -> QRectF:
        inset_x = rect.width() * scale_ratio * self._press_progress
        inset_y = rect.height() * scale_ratio * self._press_progress
        return rect.adjusted(inset_x, inset_y, -inset_x, -inset_y)

    def _update_animation_frame(self) -> None:
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        if event.button() == Qt.LeftButton:
            self._animate_press(1.0, duration=70)
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        if event.button() == Qt.LeftButton:
            self._animate_press(0.0, duration=125)
        super().mouseReleaseEvent(event)


class AcrylicButton(AnimatedButtonMixin, QPushButton):
    """玻璃质感箭头按钮。"""

    def __init__(self, text: str) -> None:
        super().__init__(text)
        self._hovered = False
        self.setCursor(Qt.PointingHandCursor)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self._init_button_animation()

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self._animate_hover(1.0)
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self._animate_hover(0.0)
        if not self.isDown():
            self._animate_press(0.0, duration=125)
        super().leaveEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        text_rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        rect = self._animated_rect(text_rect)
        path = QPainterPath()
        path.addRoundedRect(rect, 16, 16)

        top_alpha = int(_lerp(132, 188, self._hover_progress) - 30 * self._press_progress)
        bottom_alpha = int(_lerp(98, 150, self._hover_progress) - 24 * self._press_progress)

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        border_alpha = int(_lerp(178, 232, self._hover_progress))
        painter.setPen(QPen(QColor(255, 255, 255, border_alpha), 1.0))
        painter.drawPath(path)

        font = QFont(self.font())
        font.setPointSize(24)
        font.setWeight(QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.setPen(QColor("#0f172a"))
        painter.drawText(text_rect, Qt.AlignCenter, self.text())


class AcrylicTextButton(AnimatedButtonMixin, QPushButton):
    """玻璃质感文本按钮。

    COM 选择框和连接按钮共用这一套绘制逻辑，避免不同系统控件
    混用后出现高度、圆角和边框不一致的问题。
    """

    def __init__(self, text: str) -> None:
        super().__init__(text)
        self._hovered = False
        self.setCursor(Qt.PointingHandCursor)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self._init_button_animation()

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self._animate_hover(1.0)
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self._animate_hover(0.0)
        if not self.isDown():
            self._animate_press(0.0, duration=125)
        super().leaveEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        text_rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        rect = self._animated_rect(text_rect)
        path = QPainterPath()
        path.addRoundedRect(rect, 12, 12)

        active = self.isCheckable() and self.isChecked()
        if active:
            top_alpha = int(_lerp(188, 210, self._hover_progress) - 20 * self._press_progress)
            bottom_alpha = int(_lerp(146, 168, self._hover_progress) - 18 * self._press_progress)
        else:
            top_alpha = int(_lerp(154, 196, self._hover_progress) - 26 * self._press_progress)
            bottom_alpha = int(_lerp(112, 154, self._hover_progress) - 22 * self._press_progress)

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        border_alpha = int(_lerp(194, 232, max(self._hover_progress, 1.0 if active else 0.0)))
        border_color = QColor(59, 130, 246, border_alpha) if active else QColor(203, 213, 225, border_alpha)
        painter.setPen(QPen(border_color, 1.0))
        painter.drawPath(path)

        painter.setFont(self.font())
        painter.setPen(QColor("#0f172a"))
        painter.drawText(text_rect.adjusted(10, 0, -10, 0), Qt.AlignCenter, self.text())


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
        self._selection_changed = None
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

    def set_selection_changed_callback(self, callback) -> None:
        """设置选中回调，用于记住用户选择的 COM 口。"""
        self._selection_changed = callback

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
                option = AcrylicTextButton(item)
                option.setObjectName("selectOption")
                option.setFixedHeight(30)
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
        if self._selection_changed is not None:
            self._selection_changed(value)
        if self._popup is not None:
            self._popup.deleteLater()
            self._popup = None


class LedCard(AnimatedButtonMixin, QPushButton):
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
        self._current_value_ma: float | None = None
        self._decimal_places = DEFAULT_DISPLAY_DECIMALS
        self._shadow = apply_shadow(self, blur_radius=18, y_offset=4, alpha=22)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self._init_button_animation()
        self.toggled.connect(lambda _checked: self._refresh_shadow())
        self.set_current_value(None)
        self._refresh_shadow()

    def set_current_value(self, value_ma: float | None) -> None:
        """更新卡片内显示的电流值。"""
        self._current_value_ma = value_ma
        self._value_text = "--.-- mA" if value_ma is None else f"{value_ma:.{self._decimal_places}f} mA"
        self.setText(f"LED电流 {self.led_index}\n{self._value_text}")
        self.update()

    def set_decimal_places(self, decimal_places: int) -> None:
        """调整显示小数位数，并用当前值刷新文本。"""
        self._decimal_places = decimal_places
        self.set_current_value(self._current_value_ma)

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = True
        self._animate_hover(1.0)
        self._refresh_shadow()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        self._hovered = False
        self._animate_hover(0.0)
        if not self.isDown():
            self._animate_press(0.0, duration=125)
        self._refresh_shadow()
        super().leaveEvent(event)

    def _update_animation_frame(self) -> None:
        self._refresh_shadow()

    def _refresh_shadow(self) -> None:
        """根据普通、悬停、选中状态调整卡片阴影。"""
        if self.isChecked():
            blur_radius = _lerp(30, 33, self._hover_progress)
            y_offset = _lerp(9, 10, self._hover_progress)
            alpha = _lerp(46, 56, self._hover_progress)
        else:
            blur_radius = _lerp(22, 27, self._hover_progress)
            y_offset = _lerp(5, 7, self._hover_progress)
            alpha = _lerp(24, 38, self._hover_progress)
        self._shadow.setBlurRadius(blur_radius)
        self._shadow.setOffset(0, y_offset)
        if self.isChecked():
            self._shadow.setColor(QColor(59, 130, 246, int(alpha)))
        else:
            self._shadow.setColor(QColor(15, 23, 42, int(alpha)))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        text_rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        rect = self._animated_rect(text_rect, scale_ratio=0.005)
        path = QPainterPath()
        path.addRoundedRect(rect, 18, 18)

        checked = self.isChecked()
        if checked:
            top_color = QColor(219, 234, 254, int(_lerp(150, 176, self._hover_progress) - 14 * self._press_progress))
            bottom_color = QColor(219, 234, 254, int(_lerp(112, 140, self._hover_progress) - 12 * self._press_progress))
            border_color = QColor(59, 130, 246, 218)
        else:
            top_alpha = int(_lerp(128, 170, self._hover_progress) - 18 * self._press_progress)
            bottom_alpha = int(_lerp(92, 138, self._hover_progress) - 16 * self._press_progress)
            top_color = QColor(255, 255, 255, top_alpha)
            bottom_color = QColor(255, 255, 255, bottom_alpha)
            border_color = QColor(255, 255, 255, int(_lerp(178, 230, self._hover_progress)))

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
        painter.drawText(text_rect.adjusted(12, 16, -12, -text_rect.height() * 0.48), Qt.AlignHCenter | Qt.AlignTop, f"LED电流 {self.led_index}")

        painter.setFont(value_font)
        painter.drawText(text_rect.adjusted(12, text_rect.height() * 0.33, -12, -12), Qt.AlignCenter, self._value_text)


class SettingInputRow(AcrylicPanel):
    """设置页面中的单个参数行。"""

    def __init__(
        self,
        variable_name: str,
        title: str,
        description: str,
        value_text: str,
        validator: QDoubleValidator | QIntValidator,
    ) -> None:
        super().__init__(radius=14, top_alpha=128, bottom_alpha=92)
        self.setObjectName("settingRow")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(16, 13, 16, 13)
        layout.setSpacing(18)

        text_box = QWidget()
        text_box.setObjectName("settingTextBox")
        text_layout = QVBoxLayout(text_box)
        text_layout.setContentsMargins(0, 0, 0, 0)
        text_layout.setSpacing(4)

        name_label = QLabel(f"{variable_name}  {title}")
        name_label.setObjectName("settingName")
        name_label.setStyleSheet("color: #0f172a; font-weight: 600;")
        description_label = QLabel(description)
        description_label.setObjectName("settingDescription")
        description_label.setStyleSheet("color: #475569; font-size: 12px;")
        description_label.setWordWrap(True)
        text_layout.addWidget(name_label)
        text_layout.addWidget(description_label)

        self.input = QLineEdit(value_text)
        self.input.setObjectName("settingsInput")
        self.input.setValidator(validator)
        self.input.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.input.setFixedWidth(150)

        layout.addWidget(text_box, 1)
        layout.addWidget(self.input, 0, Qt.AlignRight | Qt.AlignVCenter)

    def value_text(self) -> str:
        """返回输入框文本。"""
        return self.input.text().strip()

    def set_value_text(self, value_text: str) -> None:
        """设置输入框文本并清除错误状态。"""
        self.input.setText(value_text)
        self.set_error(False)

    def set_error(self, has_error: bool) -> None:
        """切换输入框错误样式。"""
        self.input.setProperty("error", has_error)
        self.input.style().unpolish(self.input)
        self.input.style().polish(self.input)


class SettingSwitchRow(AcrylicPanel):
    """设置页面中的开关行。"""

    def __init__(
        self,
        variable_name: str,
        title: str,
        description: str,
        checked: bool,
    ) -> None:
        super().__init__(radius=14, top_alpha=128, bottom_alpha=92)
        self.setObjectName("settingRow")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(16, 13, 16, 13)
        layout.setSpacing(18)

        text_box = QWidget()
        text_layout = QVBoxLayout(text_box)
        text_layout.setContentsMargins(0, 0, 0, 0)
        text_layout.setSpacing(4)

        name_label = QLabel(f"{variable_name}  {title}")
        name_label.setObjectName("settingName")
        name_label.setStyleSheet("color: #0f172a; font-weight: 600;")
        description_label = QLabel(description)
        description_label.setObjectName("settingDescription")
        description_label.setStyleSheet("color: #475569; font-size: 12px;")
        description_label.setWordWrap(True)
        text_layout.addWidget(name_label)
        text_layout.addWidget(description_label)

        self.checkbox = QCheckBox()
        self.checkbox.setObjectName("settingsSwitch")
        self.checkbox.setChecked(checked)

        layout.addWidget(text_box, 1)
        layout.addWidget(self.checkbox, 0, Qt.AlignRight | Qt.AlignVCenter)

    def is_checked(self) -> bool:
        """返回开关状态。"""
        return self.checkbox.isChecked()

    def set_checked(self, checked: bool) -> None:
        """设置开关状态。"""
        self.checkbox.setChecked(checked)


class SettingActionRow(AcrylicPanel):
    """设置页面中的操作按钮行。"""

    def __init__(
        self,
        title: str,
        description: str,
        button_text: str,
        callback,
    ) -> None:
        super().__init__(radius=14, top_alpha=128, bottom_alpha=92)
        self.setObjectName("settingRow")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(16, 13, 16, 13)
        layout.setSpacing(18)

        text_box = QWidget()
        text_layout = QVBoxLayout(text_box)
        text_layout.setContentsMargins(0, 0, 0, 0)
        text_layout.setSpacing(4)

        name_label = QLabel(title)
        name_label.setObjectName("settingName")
        name_label.setStyleSheet("color: #0f172a; font-weight: 600;")
        description_label = QLabel(description)
        description_label.setObjectName("settingDescription")
        description_label.setStyleSheet("color: #475569; font-size: 12px;")
        description_label.setWordWrap(True)
        text_layout.addWidget(name_label)
        text_layout.addWidget(description_label)

        self.button = AcrylicTextButton(button_text)
        self.button.setFixedSize(132, 38)
        self.button.clicked.connect(callback)

        layout.addWidget(text_box, 1)
        layout.addWidget(self.button, 0, Qt.AlignRight | Qt.AlignVCenter)


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
        self.settings = load_settings()
        self.data_store = DataStore(max_points=self.settings.max_points)
        self.serial_poller: SerialPoller | None = None
        self.setWindowTitle("LED亮度检测上位机")
        self.setWindowIcon(QIcon(str(APP_ICON_PATH)))
        self.resize(1080, 720)
        self._build_ui()
        self._apply_display_decimals()
        apply_native_title_bar_style(self)
        self._select_led(1)
        self._refresh_ports()
        QTimer.singleShot(0, self._auto_connect_last_port)

    def showEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        super().showEvent(event)
        apply_native_title_bar_style(self)

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
        self.port_combo.set_selection_changed_callback(self._handle_port_selected)

        self.connect_button = AcrylicTextButton("连接")
        self.connect_button.setObjectName("connectButton")
        self.connect_button.setFixedSize(92, 38)
        self.connect_button.clicked.connect(self._toggle_connection)

        self.status_label = QLabel("状态：未连接")
        self.status_label.setObjectName("statusLabel")

        self.settings_button = AcrylicTextButton("设置")
        self.settings_button.setObjectName("settingsButton")
        self.settings_button.setFixedSize(92, 38)
        self.settings_button.clicked.connect(self._open_settings)

        toolbar.addWidget(self.port_combo)
        toolbar.addWidget(self.connect_button)
        toolbar.addWidget(self.status_label)
        toolbar.addItem(QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
        toolbar.addWidget(self.settings_button)

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
        self.plot_widget.setBackground(None)
        self.plot_widget.setStyleSheet("background: transparent; border: none;")
        self.plot_widget.viewport().setStyleSheet("background: transparent;")
        self.plot_widget.showGrid(x=self.settings.curve_show_grid, y=self.settings.curve_show_grid, alpha=0.08)
        self.plot_widget.setLabel("left", "电流值", units="mA")
        self.plot_widget.setLabel("bottom", "时间")
        self.plot_widget.setMouseEnabled(x=False, y=False)
        self.plot_widget.getPlotItem().setMenuEnabled(False)
        self.plot_widget.getPlotItem().hideButtons()
        self.plot_widget.getAxis("left").setPen(pg.mkPen("#0f172a"))
        self.plot_widget.getAxis("bottom").setPen(pg.mkPen("#0f172a"))
        self.plot_widget.getAxis("left").setTextPen(pg.mkPen("#0f172a"))
        self.plot_widget.getAxis("bottom").setTextPen(pg.mkPen("#0f172a"))
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

        curve_layout.addWidget(self.curve_title)
        curve_layout.addWidget(self.plot_widget, 1)

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

        self.monitor_view = content
        self.settings_view = self._build_settings_view()
        self.page_stack = QStackedWidget()
        self.page_stack.setObjectName("pageStack")
        self.page_stack.addWidget(self.monitor_view)
        self.page_stack.addWidget(self.settings_view)
        root_layout.addWidget(self.page_stack, 1)
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
            #pageStack {
                background: transparent;
            }
            #settingsPageStack {
                background: transparent;
            }
            #settingsScroll {
                background: transparent;
                border: none;
            }
            #settingsScroll > QWidget > QWidget {
                background: transparent;
            }
            #settingsScroll QScrollBar:vertical {
                width: 8px;
                margin: 4px 0 4px 0;
                border: none;
                border-radius: 4px;
                background: rgba(226, 232, 240, 0.28);
            }
            #settingsScroll QScrollBar::handle:vertical {
                min-height: 44px;
                border-radius: 4px;
                background: rgba(96, 165, 250, 0.44);
            }
            #settingsScroll QScrollBar::handle:vertical:hover {
                background: rgba(59, 130, 246, 0.58);
            }
            #settingsScroll QScrollBar::add-line:vertical,
            #settingsScroll QScrollBar::sub-line:vertical {
                height: 0;
                border: none;
                background: transparent;
            }
            #settingsScroll QScrollBar::add-page:vertical,
            #settingsScroll QScrollBar::sub-page:vertical {
                background: transparent;
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
            #settingsSideTitle {
                color: #0f172a;
                font-size: 26px;
                font-weight: 600;
            }
            #settingsSideDescription {
                color: #64748b;
                font-size: 13px;
                line-height: 1.5;
            }
            #settingsSectionTitle {
                color: #0f172a;
                font-size: 18px;
                font-weight: 600;
            }
            #settingsName {
                color: #0f172a;
                font-size: 14px;
                font-weight: 600;
            }
            #settingsDescription {
                color: #64748b;
                font-size: 12px;
            }
            #settingsInput {
                min-height: 32px;
                border: 1px solid rgba(148, 163, 184, 0.62);
                border-radius: 8px;
                padding: 0 10px;
                background: rgba(255, 255, 255, 0.78);
                color: #0f172a;
                selection-background-color: #bfdbfe;
            }
            #settingsInput:focus {
                border-color: #60a5fa;
                background: rgba(255, 255, 255, 0.94);
            }
            #settingsInput[error="true"] {
                border-color: #ef4444;
                background: rgba(254, 242, 242, 0.88);
            }
            #settingsSwitch {
                background: transparent;
            }
            #settingsSwitch::indicator {
                width: 38px;
                height: 22px;
                border-radius: 11px;
                border: 1px solid rgba(148, 163, 184, 0.86);
                background: rgba(255, 255, 255, 0.72);
            }
            #settingsSwitch::indicator:checked {
                border-color: rgba(59, 130, 246, 0.92);
                background: rgba(96, 165, 250, 0.78);
            }
            #settingsHint {
                color: #64748b;
                font-size: 12px;
            }
            """
        )

    def _build_settings_view(self) -> QWidget:
        """构建嵌入主窗口的设置页面。"""
        settings_panel = AcrylicPanel(radius=26, top_alpha=142, bottom_alpha=104)
        settings_panel.setObjectName("settingsPanel")

        layout = QHBoxLayout(settings_panel)
        layout.setContentsMargins(24, 22, 24, 22)
        layout.setSpacing(28)

        side = QWidget()
        side.setObjectName("settingsSide")
        side.setFixedWidth(250)
        side_layout = QVBoxLayout(side)
        side_layout.setContentsMargins(0, 0, 0, 0)
        side_layout.setSpacing(12)

        title = QLabel("设置")
        title.setObjectName("settingsSideTitle")
        self.common_settings_button = AcrylicTextButton("常用设置")
        self.common_settings_button.setFixedSize(178, 38)
        self.common_settings_button.setCheckable(True)
        self.common_settings_button.clicked.connect(lambda: self._set_settings_section(0))
        self.parameter_settings_button = AcrylicTextButton("采集与曲线设置")
        self.parameter_settings_button.setFixedSize(178, 38)
        self.parameter_settings_button.setCheckable(True)
        self.parameter_settings_button.clicked.connect(lambda: self._set_settings_section(1))

        side_layout.addWidget(title)
        side_layout.addSpacing(12)
        side_layout.addWidget(self.common_settings_button)
        side_layout.addWidget(self.parameter_settings_button)
        side_layout.addStretch(1)

        settings_content = QWidget()
        settings_content.setObjectName("settingsContent")
        content_layout = QVBoxLayout(settings_content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(12)

        self.settings_page_stack = QStackedWidget()
        self.settings_page_stack.setObjectName("settingsPageStack")
        self.settings_page_stack.addWidget(self._wrap_settings_scroll(self._build_common_settings_page()))
        self.settings_page_stack.addWidget(self._wrap_settings_scroll(self._build_parameter_settings_page()))
        content_layout.addWidget(self.settings_page_stack, 1)

        self.settings_hint = QLabel("")
        self.settings_hint.setObjectName("settingsHint")
        self.settings_hint.setWordWrap(True)
        content_layout.addWidget(self.settings_hint)

        button_layout = QHBoxLayout()
        button_layout.setContentsMargins(0, 8, 0, 0)
        button_layout.addStretch(1)

        self.restore_defaults_button = AcrylicTextButton("恢复默认")
        self.restore_defaults_button.setFixedSize(104, 38)
        self.restore_defaults_button.clicked.connect(self._restore_default_settings_inputs)

        self.discard_settings_button = AcrylicTextButton("放弃更改")
        self.discard_settings_button.setFixedSize(104, 38)
        self.discard_settings_button.clicked.connect(self._discard_settings_changes)

        self.save_settings_button = AcrylicTextButton("保存")
        self.save_settings_button.setFixedSize(92, 38)
        self.save_settings_button.clicked.connect(self._save_settings_from_page)

        button_layout.addWidget(self.restore_defaults_button)
        button_layout.addWidget(self.discard_settings_button)
        button_layout.addWidget(self.save_settings_button)
        content_layout.addLayout(button_layout)

        layout.addWidget(side)
        layout.addWidget(settings_content, 1)
        self._set_settings_section(0)
        return settings_panel

    def _wrap_settings_scroll(self, page: QWidget) -> QScrollArea:
        """为设置页添加透明滚动容器。"""
        scroll_area = QScrollArea()
        scroll_area.setObjectName("settingsScroll")
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.viewport().setStyleSheet("background: transparent;")
        scroll_area.setWidget(page)
        return scroll_area

    def _build_parameter_settings_page(self) -> QWidget:
        """构建采集与曲线参数设置页。"""
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 14, 0)
        layout.setSpacing(12)

        section_title = QLabel("采集与曲线设置")
        section_title.setObjectName("settingsSectionTitle")
        layout.addWidget(section_title)

        read_timeout_validator = QDoubleValidator(READ_TIMEOUT_RANGE[0], READ_TIMEOUT_RANGE[1], 2, self)
        read_timeout_validator.setNotation(QDoubleValidator.Notation.StandardNotation)
        round_delay_validator = QDoubleValidator(ROUND_DELAY_RANGE[0], ROUND_DELAY_RANGE[1], 2, self)
        round_delay_validator.setNotation(QDoubleValidator.Notation.StandardNotation)
        curve_seconds_validator = QIntValidator(
            int(CURVE_VISIBLE_SECONDS_RANGE[0]),
            int(CURVE_VISIBLE_SECONDS_RANGE[1]),
            self,
        )
        max_points_validator = QIntValidator(MAX_POINTS_RANGE[0], MAX_POINTS_RANGE[1], self)

        self.read_timeout_row = SettingInputRow(
            "READ_TIMEOUT_SECONDS",
            "读取超时",
            "串口读写返回帧的最大等待时间，单位秒。过小可能误报超时，过大异常响应会变慢。",
            self._format_float(self.settings.read_timeout_seconds),
            read_timeout_validator,
        )
        self.round_delay_row = SettingInputRow(
            "ROUND_DELAY_SECONDS",
            "轮询间隔",
            "完成 1-16 通道一轮读取后的暂停时间，单位秒。数值越小刷新越快。",
            self._format_float(self.settings.round_delay_seconds),
            round_delay_validator,
        )
        self.curve_seconds_row = SettingInputRow(
            "CURVE_VISIBLE_SECONDS",
            "曲线显示时长",
            "右侧曲线固定显示最近多少秒的数据，只影响显示窗口，不改变采集。",
            self._format_int(self.settings.curve_visible_seconds),
            curve_seconds_validator,
        )
        self.max_points_row = SettingInputRow(
            "DEFAULT_MAX_POINTS",
            "历史最大点数",
            "每个通道/LED 在内存中保留的历史点数量，用于控制长时间运行的内存占用。",
            self._format_int(self.settings.max_points),
            max_points_validator,
        )

        for row in (
            self.read_timeout_row,
            self.round_delay_row,
            self.curve_seconds_row,
            self.max_points_row,
        ):
            self._add_settings_row(layout, row)

        layout.addStretch(1)
        return page

    def _build_common_settings_page(self) -> QWidget:
        """构建常用设置页。"""
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 14, 0)
        layout.setSpacing(12)

        section_title = QLabel("常用设置")
        section_title.setObjectName("settingsSectionTitle")
        layout.addWidget(section_title)

        decimals_validator = QIntValidator(DISPLAY_DECIMALS_RANGE[0], DISPLAY_DECIMALS_RANGE[1], self)
        line_width_validator = QDoubleValidator(CURVE_LINE_WIDTH_RANGE[0], CURVE_LINE_WIDTH_RANGE[1], 1, self)
        line_width_validator.setNotation(QDoubleValidator.Notation.StandardNotation)
        y_axis_validator = QDoubleValidator(Y_AXIS_RANGE[0], Y_AXIS_RANGE[1], 3, self)
        y_axis_validator.setNotation(QDoubleValidator.Notation.StandardNotation)

        self.remember_port_row = SettingSwitchRow(
            "REMEMBER_LAST_PORT",
            "记住上次选择的 COM 口",
            "开启后会保存最近选择或连接的串口，下次启动时自动选中。",
            self.settings.remember_last_port,
        )
        self.auto_connect_row = SettingSwitchRow(
            "AUTO_CONNECT_LAST_PORT",
            "启动后自动连接上次串口",
            "开启后程序启动时会尝试连接已记住的串口。默认关闭，避免误连接。",
            self.settings.auto_connect_last_port,
        )
        self.display_decimals_row = SettingInputRow(
            "DISPLAY_DECIMALS",
            "显示小数位数",
            "LED 电流卡片显示的小数位数量。范围 0 到 6。",
            self._format_int(self.settings.display_decimals),
            decimals_validator,
        )
        self.curve_grid_row = SettingSwitchRow(
            "CURVE_SHOW_GRID",
            "显示曲线网格",
            "控制右侧曲线区域是否显示浅色网格线。",
            self.settings.curve_show_grid,
        )
        self.curve_line_width_row = SettingInputRow(
            "CURVE_LINE_WIDTH",
            "曲线线宽",
            "实时曲线线条宽度。范围 1.0 到 6.0，实际绘制时保留整数线宽。",
            self._format_float(self.settings.curve_line_width),
            line_width_validator,
        )
        self.y_axis_auto_row = SettingSwitchRow(
            "Y_AXIS_AUTO_SCALE",
            "Y 轴自动缩放",
            "开启时曲线纵轴随数据自动调整；关闭时使用下面的固定范围。",
            self.settings.y_axis_auto_scale,
        )
        self.y_axis_min_row = SettingInputRow(
            "Y_AXIS_MIN",
            "Y 轴固定下限",
            "关闭自动缩放后使用的纵轴最小值。",
            self._format_float(self.settings.y_axis_min),
            y_axis_validator,
        )
        self.y_axis_max_row = SettingInputRow(
            "Y_AXIS_MAX",
            "Y 轴固定上限",
            "关闭自动缩放后使用的纵轴最大值。",
            self._format_float(self.settings.y_axis_max),
            y_axis_validator,
        )
        self.show_serial_errors_row = SettingSwitchRow(
            "SHOW_SERIAL_ERRORS",
            "串口读取失败提示到状态栏",
            "开启后读取失败会刷新顶部状态栏；关闭后减少频繁错误提示干扰。",
            self.settings.show_serial_errors,
        )
        clear_history_row = SettingActionRow(
            "清空当前曲线历史",
            "只清空当前通道、当前 LED 的内存曲线数据，不影响最新读数和其他通道。",
            "清空当前曲线",
            self._clear_current_curve_history,
        )

        for row in (
            self.remember_port_row,
            self.auto_connect_row,
            self.display_decimals_row,
            self.curve_grid_row,
            self.curve_line_width_row,
            self.y_axis_auto_row,
            self.y_axis_min_row,
            self.y_axis_max_row,
            self.show_serial_errors_row,
            clear_history_row,
        ):
            self._add_settings_row(layout, row)

        layout.addStretch(1)
        return page

    def _add_settings_row(self, layout: QVBoxLayout, row: QWidget) -> None:
        """把设置行添加到布局并套用轻量阴影。"""
        apply_shadow(row, blur_radius=18, y_offset=5, alpha=18)
        layout.addWidget(row)

    def _set_settings_section(self, index: int) -> None:
        """切换设置页右侧具体分类。"""
        self.settings_page_stack.setCurrentIndex(index)
        self.common_settings_button.setChecked(index == 0)
        self.parameter_settings_button.setChecked(index == 1)
        if index == 0:
            self.settings_hint.setText("常用设置保存后立即应用；自动连接会在下一次启动时生效。")
        else:
            self.settings_hint.setText("读取超时和轮询间隔会在下一次连接串口时生效；曲线显示时长和历史点数保存后立即生效。")

    def _format_float(self, value: float) -> str:
        """格式化设置页面中的浮点数，避免不必要的尾随零。"""
        return f"{value:.2f}".rstrip("0").rstrip(".")

    def _format_int(self, value: float | int) -> str:
        """格式化设置页面中的整数。"""
        return str(int(value))

    def _populate_settings_inputs(self, settings: AppSettings) -> None:
        """用当前设置刷新设置页输入框。"""
        self.read_timeout_row.set_value_text(self._format_float(settings.read_timeout_seconds))
        self.round_delay_row.set_value_text(self._format_float(settings.round_delay_seconds))
        self.curve_seconds_row.set_value_text(self._format_int(settings.curve_visible_seconds))
        self.max_points_row.set_value_text(self._format_int(settings.max_points))
        self.remember_port_row.set_checked(settings.remember_last_port)
        self.auto_connect_row.set_checked(settings.auto_connect_last_port)
        self.display_decimals_row.set_value_text(self._format_int(settings.display_decimals))
        self.curve_grid_row.set_checked(settings.curve_show_grid)
        self.curve_line_width_row.set_value_text(self._format_float(settings.curve_line_width))
        self.y_axis_auto_row.set_checked(settings.y_axis_auto_scale)
        self.y_axis_min_row.set_value_text(self._format_float(settings.y_axis_min))
        self.y_axis_max_row.set_value_text(self._format_float(settings.y_axis_max))
        self.show_serial_errors_row.set_checked(settings.show_serial_errors)

    def _show_monitor_page(self, status_text: str | None = None) -> None:
        """切回监测页面。"""
        self.page_stack.setCurrentWidget(self.monitor_view)
        self.settings_button.setText("设置")
        if status_text is not None:
            self.status_label.setText(status_text)

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
        if (
            self.settings.remember_last_port
            and self.port_combo.currentText() not in self.available_ports
            and self.settings.last_port_name in self.available_ports
        ):
            self.port_combo.setText(self.settings.last_port_name)
        if self.serial_poller is None:
            if self.available_ports:
                self.status_label.setText(f"状态：未连接，发现 {len(self.available_ports)} 个COM口")
            else:
                self.status_label.setText("状态：未连接，未发现COM口")

    def _handle_port_selected(self, port_name: str) -> None:
        """记住用户手动选择的 COM 口。"""
        if not self.settings.remember_last_port:
            return
        self.settings = replace(self.settings, last_port_name=port_name)
        save_settings(self.settings)

    def _auto_connect_last_port(self) -> None:
        """启动后按设置自动连接上次串口。"""
        if (
            not self.settings.auto_connect_last_port
            or not self.settings.remember_last_port
            or self.serial_poller is not None
            or self.settings.last_port_name not in self.available_ports
        ):
            return
        self.port_combo.setText(self.settings.last_port_name)
        self._toggle_connection()

    def _open_settings(self) -> None:
        """在主窗口内切换到设置页面。"""
        if self.page_stack.currentWidget() == self.settings_view:
            self._discard_settings_changes()
            return
        self._populate_settings_inputs(self.settings)
        self.page_stack.setCurrentWidget(self.settings_view)
        self.settings_button.setText("返回")
        self.status_label.setText("状态：正在修改设置")

    def _restore_default_settings_inputs(self) -> None:
        """把设置页输入框恢复为默认值，点击保存后才会生效。"""
        self._populate_settings_inputs(
            AppSettings(
                read_timeout_seconds=DEFAULT_READ_TIMEOUT_SECONDS,
                round_delay_seconds=DEFAULT_ROUND_DELAY_SECONDS,
                curve_visible_seconds=DEFAULT_CURVE_VISIBLE_SECONDS,
                max_points=DEFAULT_MAX_POINTS,
                display_decimals=DEFAULT_DISPLAY_DECIMALS,
                curve_show_grid=DEFAULT_CURVE_SHOW_GRID,
                curve_line_width=DEFAULT_CURVE_LINE_WIDTH,
                y_axis_auto_scale=DEFAULT_Y_AXIS_AUTO_SCALE,
                y_axis_min=DEFAULT_Y_AXIS_MIN,
                y_axis_max=DEFAULT_Y_AXIS_MAX,
                show_serial_errors=DEFAULT_SHOW_SERIAL_ERRORS,
            )
        )
        self.status_label.setText("状态：已填入默认设置，点击保存后生效")

    def _discard_settings_changes(self) -> None:
        """放弃未保存的输入并回到监测页面。"""
        self._populate_settings_inputs(self.settings)
        self._show_monitor_page("状态：已放弃设置更改")

    def _save_settings_from_page(self) -> None:
        """校验设置页输入并保存。"""
        new_settings = self._settings_from_inputs()
        if new_settings is None:
            return

        port_name = self.port_combo.currentText()
        if new_settings.remember_last_port and port_name in self.available_ports:
            new_settings = replace(new_settings, last_port_name=port_name)
        self.settings = new_settings
        save_settings(self.settings)
        self._apply_runtime_settings()
        self._update_curve()
        if self.serial_poller is None:
            self._show_monitor_page("状态：设置已保存")
        else:
            self._show_monitor_page("状态：设置已保存，串口参数下次连接生效")

    def _settings_from_inputs(self) -> AppSettings | None:
        """从设置页输入框读取并校验参数。"""
        for row in (
            self.read_timeout_row,
            self.round_delay_row,
            self.curve_seconds_row,
            self.max_points_row,
            self.display_decimals_row,
            self.curve_line_width_row,
            self.y_axis_min_row,
            self.y_axis_max_row,
        ):
            row.set_error(False)

        read_timeout = self._parse_float_setting(
            self.read_timeout_row,
            READ_TIMEOUT_RANGE,
            "READ_TIMEOUT_SECONDS",
        )
        round_delay = self._parse_float_setting(
            self.round_delay_row,
            ROUND_DELAY_RANGE,
            "ROUND_DELAY_SECONDS",
        )
        curve_seconds = self._parse_int_setting(
            self.curve_seconds_row,
            (int(CURVE_VISIBLE_SECONDS_RANGE[0]), int(CURVE_VISIBLE_SECONDS_RANGE[1])),
            "CURVE_VISIBLE_SECONDS",
        )
        max_points = self._parse_int_setting(
            self.max_points_row,
            MAX_POINTS_RANGE,
            "DEFAULT_MAX_POINTS",
        )
        display_decimals = self._parse_int_setting(
            self.display_decimals_row,
            DISPLAY_DECIMALS_RANGE,
            "DISPLAY_DECIMALS",
        )
        curve_line_width = self._parse_float_setting(
            self.curve_line_width_row,
            CURVE_LINE_WIDTH_RANGE,
            "CURVE_LINE_WIDTH",
        )
        y_axis_min = self._parse_float_setting(
            self.y_axis_min_row,
            Y_AXIS_RANGE,
            "Y_AXIS_MIN",
        )
        y_axis_max = self._parse_float_setting(
            self.y_axis_max_row,
            Y_AXIS_RANGE,
            "Y_AXIS_MAX",
        )

        if None in (
            read_timeout,
            round_delay,
            curve_seconds,
            max_points,
            display_decimals,
            curve_line_width,
            y_axis_min,
            y_axis_max,
        ):
            return None
        if y_axis_min >= y_axis_max:
            self.y_axis_min_row.set_error(True)
            self.y_axis_max_row.set_error(True)
            self.status_label.setText("状态：Y_AXIS_MIN 必须小于 Y_AXIS_MAX")
            self.y_axis_min_row.input.setFocus()
            return None
        return AppSettings(
            read_timeout_seconds=read_timeout,
            round_delay_seconds=round_delay,
            curve_visible_seconds=float(curve_seconds),
            max_points=max_points,
            remember_last_port=self.remember_port_row.is_checked(),
            last_port_name=self.settings.last_port_name,
            auto_connect_last_port=self.auto_connect_row.is_checked(),
            display_decimals=display_decimals,
            curve_show_grid=self.curve_grid_row.is_checked(),
            curve_line_width=curve_line_width,
            y_axis_auto_scale=self.y_axis_auto_row.is_checked(),
            y_axis_min=y_axis_min,
            y_axis_max=y_axis_max,
            show_serial_errors=self.show_serial_errors_row.is_checked(),
        )

    def _parse_float_setting(
        self,
        row: SettingInputRow,
        value_range: tuple[float, float],
        variable_name: str,
    ) -> float | None:
        """校验浮点参数输入。"""
        try:
            value = float(row.value_text())
        except ValueError:
            row.set_error(True)
            self.status_label.setText(f"状态：{variable_name} 请输入数字")
            row.input.setFocus()
            return None
        if not value_range[0] <= value <= value_range[1]:
            row.set_error(True)
            self.status_label.setText(f"状态：{variable_name} 范围为 {value_range[0]} 到 {value_range[1]}")
            row.input.setFocus()
            return None
        return value

    def _parse_int_setting(
        self,
        row: SettingInputRow,
        value_range: tuple[int, int],
        variable_name: str,
    ) -> int | None:
        """校验整数参数输入。"""
        try:
            value = int(row.value_text())
        except ValueError:
            row.set_error(True)
            self.status_label.setText(f"状态：{variable_name} 请输入整数")
            row.input.setFocus()
            return None
        if not value_range[0] <= value <= value_range[1]:
            row.set_error(True)
            self.status_label.setText(f"状态：{variable_name} 范围为 {value_range[0]} 到 {value_range[1]}")
            row.input.setFocus()
            return None
        return value

    def _apply_runtime_settings(self) -> None:
        """把已保存设置应用到当前运行中的界面状态。"""
        self.data_store.set_max_points(self.settings.max_points)
        self._apply_display_decimals()
        self._apply_curve_settings()

    def _apply_display_decimals(self) -> None:
        """应用 LED 卡片显示小数位设置。"""
        for card in self.led_cards:
            card.set_decimal_places(self.settings.display_decimals)

    def _apply_curve_settings(self) -> None:
        """应用曲线显示设置。"""
        self.plot_widget.showGrid(
            x=self.settings.curve_show_grid,
            y=self.settings.curve_show_grid,
            alpha=0.08,
        )
        self.plot_curve.setPen(self._curve_pen())
        x_values, y_values = self.plot_curve.getData()
        if x_values is None or y_values is None:
            self.plot_curve.setData([], [], antialias=False, stepMode=CURVE_STEP_MODE)
            self._apply_curve_y_range([])
        else:
            self.plot_curve.setData(x_values, y_values, antialias=False, stepMode=CURVE_STEP_MODE)
            self._apply_curve_y_range(list(y_values))

    def _curve_pen(self):
        """返回整数像素线宽的曲线画笔，减少亚像素抗锯齿伪影。"""
        return pg.mkPen("#2563eb", width=max(1, round(self.settings.curve_line_width)))

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

    def _clear_current_curve_history(self) -> None:
        """清空当前通道、当前 LED 的内存曲线历史。"""
        self.data_store.clear_history(self.current_channel, self.selected_led)
        self._update_curve()
        self.status_label.setText(f"状态：已清空通道 {self.current_channel} / LED {self.selected_led} 曲线历史")

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

        if self.settings.remember_last_port:
            self.settings = replace(self.settings, last_port_name=port_name)
            save_settings(self.settings)

        self.serial_poller = SerialPoller(
            port_name,
            self,
            read_timeout_seconds=self.settings.read_timeout_seconds,
            round_delay_seconds=self.settings.round_delay_seconds,
        )
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
        if not self.settings.show_serial_errors:
            return
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

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        """关闭窗口前停止后台串口线程，避免进程残留。"""
        if self.serial_poller is not None:
            self.serial_poller.stop()
            self.serial_poller.wait(1000)
        super().closeEvent(event)


def run() -> int:
    """创建 Qt 应用并显示主窗口。"""
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon(str(APP_ICON_PATH)))
    window = MainWindow()
    window.show()
    return app.exec()
