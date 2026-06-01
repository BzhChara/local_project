from __future__ import annotations

from PySide6.QtCore import QEasingCurve, QEvent, QPropertyAnimation, QRectF, QSize, Qt, Property
from PySide6.QtGui import (
    QColor,
    QDoubleValidator,
    QFont,
    QIntValidator,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
)
from PySide6.QtWidgets import (
    QCheckBox,
    QFrame,
    QGraphicsDropShadowEffect,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from led_debugger.settings import DEFAULT_DISPLAY_DECIMALS


def apply_shadow(widget: QWidget, blur_radius: float, y_offset: float, alpha: int) -> QGraphicsDropShadowEffect:
    """给控件添加轻量阴影，用来模拟 Fluent/Acrylic 的悬浮层次。"""
    shadow = QGraphicsDropShadowEffect(widget)
    shadow.setBlurRadius(blur_radius)
    shadow.setOffset(0, y_offset)
    shadow.setColor(QColor(15, 23, 42, alpha))
    widget.setGraphicsEffect(shadow)
    return shadow


class ElidedLabel(QLabel):
    """宽度不足时用 ... 省略显示，完整内容保留在 tooltip。"""

    def __init__(self, text: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._full_text = ""
        self.setMinimumWidth(0)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.setText(text)

    def setText(self, text: str) -> None:  # noqa: N802 - Qt 接口沿用 Qt 命名。
        self._full_text = text
        self.setToolTip(text)
        self._refresh_elided_text()

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        super().resizeEvent(event)
        self._refresh_elided_text()

    def changeEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        super().changeEvent(event)
        if event.type() in (QEvent.Type.FontChange, QEvent.Type.StyleChange):
            self._refresh_elided_text()

    def sizeHint(self) -> QSize:  # noqa: N802 - Qt 接口沿用 Qt 命名。
        hint = super().sizeHint()
        return QSize(min(hint.width(), 260), hint.height())

    def minimumSizeHint(self) -> QSize:  # noqa: N802 - Qt 接口沿用 Qt 命名。
        hint = super().minimumSizeHint()
        return QSize(0, hint.height())

    def _refresh_elided_text(self) -> None:
        available_width = self.contentsRect().width()
        QLabel.setText(self, _elide_right(self._full_text, available_width, self.fontMetrics()))


def _elide_right(text: str, available_width: int, metrics) -> str:
    suffix = "..."
    if available_width <= 0:
        return ""
    if metrics.horizontalAdvance(text) <= available_width:
        return text

    suffix_width = metrics.horizontalAdvance(suffix)
    if suffix_width >= available_width:
        return suffix if suffix_width <= available_width else ""

    low = 0
    high = len(text)
    while low < high:
        middle = (low + high + 1) // 2
        if metrics.horizontalAdvance(text[:middle]) + suffix_width <= available_width:
            low = middle
        else:
            high = middle - 1
    return text[:low] + suffix


class AcrylicPanel(QFrame):
    """可复用的伪亚克力面板。"""

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
    """玻璃质感文本按钮。"""

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
    """自定义 COM 选择按钮。"""

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
    """单个 LED 电流显示卡片。"""

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

        self.name_label = QLabel(f"{variable_name}  {title}")
        self.name_label.setObjectName("settingName")
        self.description_label = QLabel(description)
        self.description_label.setObjectName("settingDescription")
        self.description_label.setWordWrap(True)
        text_layout.addWidget(self.name_label)
        text_layout.addWidget(self.description_label)

        self.input = QLineEdit(value_text)
        self.input.setObjectName("settingsInput")
        self.input.setValidator(validator)
        self.input.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.input.setFixedWidth(150)

        layout.addWidget(text_box, 1)
        layout.addWidget(self.input, 0, Qt.AlignRight | Qt.AlignVCenter)
        self.set_row_enabled(True)

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

    def set_row_enabled(self, enabled: bool) -> None:
        """切换整行手动输入是否可用。"""
        self.input.setEnabled(enabled)
        self.name_label.setStyleSheet(
            f"color: {'#0f172a' if enabled else '#94a3b8'}; font-weight: 600;"
        )
        self.description_label.setStyleSheet(
            f"color: {'#475569' if enabled else '#94a3b8'}; font-size: 12px;"
        )
        self.set_error(False)


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

        self.name_label = QLabel(f"{variable_name}  {title}")
        self.name_label.setObjectName("settingName")
        self.description_label = QLabel(description)
        self.description_label.setObjectName("settingDescription")
        self.description_label.setWordWrap(True)
        text_layout.addWidget(self.name_label)
        text_layout.addWidget(self.description_label)

        self.checkbox = ToggleSwitch()
        self.checkbox.setObjectName("settingsSwitch")
        self.checkbox.set_instant_checked(checked)

        layout.addWidget(text_box, 1)
        layout.addWidget(self.checkbox, 0, Qt.AlignRight | Qt.AlignVCenter)
        self.set_row_enabled(True)

    def is_checked(self) -> bool:
        """返回开关状态。"""
        return self.checkbox.isChecked()

    def set_checked(self, checked: bool) -> None:
        """设置开关状态。"""
        self.checkbox.set_instant_checked(checked)

    def set_row_enabled(self, enabled: bool) -> None:
        """切换整行开关是否可用。"""
        self.checkbox.setEnabled(enabled)
        self.name_label.setStyleSheet(
            f"color: {'#0f172a' if enabled else '#94a3b8'}; font-weight: 600;"
        )
        self.description_label.setStyleSheet(
            f"color: {'#475569' if enabled else '#94a3b8'}; font-size: 12px;"
        )


class ToggleSwitch(QCheckBox):
    """设置页使用的轻量滑动开关。"""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._thumb_position = 1.0 if self.isChecked() else 0.0
        self.setFixedSize(44, 24)
        self.setCursor(Qt.PointingHandCursor)
        self.setFocusPolicy(Qt.StrongFocus)

        self._animation = QPropertyAnimation(self, b"thumbPosition", self)
        self._animation.setDuration(130)
        self._animation.setEasingCurve(QEasingCurve.Type.OutCubic)
        self.toggled.connect(self._start_toggle_animation)

    def sizeHint(self):  # noqa: N802 - Qt 命名
        return self.size()

    def set_instant_checked(self, checked: bool) -> None:
        """程序刷新设置值时立即同步滑块位置，不播放动画。"""
        previous_blocked = self.blockSignals(True)
        self.setChecked(checked)
        self.blockSignals(previous_blocked)
        self._animation.stop()
        self._thumb_position = 1.0 if checked else 0.0
        self.update()

    def hitButton(self, pos) -> bool:  # noqa: N802 - Qt 命名
        return self.rect().contains(pos)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        track_rect = QRectF(2, 3, self.width() - 4, self.height() - 6)
        if self.isChecked():
            track_color = QColor("#334155")
            border_color = QColor("#475569")
            thumb_color = QColor("#f8fafc")
        else:
            track_color = QColor("#d1d5db")
            border_color = QColor("#cbd5e1")
            thumb_color = QColor("#0f172a")

        if not self.isEnabled():
            track_color.setAlpha(118)
            border_color.setAlpha(118)
            thumb_color.setAlpha(128)

        painter.setPen(QPen(border_color, 1.0))
        painter.setBrush(track_color)
        painter.drawRoundedRect(track_rect, track_rect.height() / 2, track_rect.height() / 2)

        thumb_size = 14.0
        left = track_rect.left() + 3.0
        right = track_rect.right() - thumb_size - 3.0
        thumb_x = left + (right - left) * self._thumb_position
        thumb_rect = QRectF(thumb_x, track_rect.center().y() - thumb_size / 2, thumb_size, thumb_size)
        painter.setPen(Qt.NoPen)
        painter.setBrush(thumb_color)
        painter.drawEllipse(thumb_rect)

    def _start_toggle_animation(self, checked: bool) -> None:
        self._animation.stop()
        self._animation.setStartValue(self._thumb_position)
        self._animation.setEndValue(1.0 if checked else 0.0)
        self._animation.start()

    def _get_thumb_position(self) -> float:
        return self._thumb_position

    def _set_thumb_position(self, value: float) -> None:
        self._thumb_position = value
        self.update()

    thumbPosition = Property(float, _get_thumb_position, _set_thumb_position)


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
