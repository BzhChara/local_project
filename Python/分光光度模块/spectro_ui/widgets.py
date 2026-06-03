from __future__ import annotations

from PySide6.QtCore import QEasingCurve, QEvent, QPropertyAnimation, QRectF, QSize, Qt, Property
from PySide6.QtGui import QColor, QFont, QLinearGradient, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QFrame, QGraphicsDropShadowEffect, QLabel, QPushButton, QSizePolicy, QVBoxLayout, QWidget


DOT_COLORS = {
    "idle": QColor("#94a3b8"),
    "ready": QColor("#22c55e"),
    "detecting": QColor("#f59e0b"),
    "waiting_air": QColor("#3b82f6"),
}


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
        if not self.isEnabled():
            top_alpha = 96
            bottom_alpha = 72
            border_color = QColor(203, 213, 225, 110)
            text_color = QColor("#94a3b8")
        elif active:
            top_alpha = int(_lerp(188, 210, self._hover_progress) - 20 * self._press_progress)
            bottom_alpha = int(_lerp(146, 168, self._hover_progress) - 18 * self._press_progress)
            border_alpha = int(_lerp(194, 232, max(self._hover_progress, 1.0)))
            border_color = QColor(59, 130, 246, border_alpha)
            text_color = QColor("#0f172a")
        else:
            top_alpha = int(_lerp(154, 196, self._hover_progress) - 26 * self._press_progress)
            bottom_alpha = int(_lerp(112, 154, self._hover_progress) - 22 * self._press_progress)
            border_alpha = int(_lerp(194, 232, self._hover_progress))
            border_color = QColor(203, 213, 225, border_alpha)
            text_color = QColor("#0f172a")

        fill = QLinearGradient(rect.topLeft(), rect.bottomLeft())
        fill.setColorAt(0.0, QColor(255, 255, 255, top_alpha))
        fill.setColorAt(1.0, QColor(255, 255, 255, bottom_alpha))
        painter.setPen(Qt.NoPen)
        painter.setBrush(fill)
        painter.drawPath(path)

        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(border_color, 1.0))
        painter.drawPath(path)

        painter.setFont(self.font())
        painter.setPen(text_color)
        painter.drawText(text_rect.adjusted(10, 0, -10, 0), Qt.AlignCenter, self.text())


class AcrylicSelectButton(AcrylicTextButton):
    """自定义选择按钮，沿用参考上位机的内部浮层样式。"""

    def __init__(self, placeholder: str) -> None:
        super().__init__(placeholder)
        self._placeholder = placeholder
        self._items: list[tuple[str, object]] = []
        self._current_index = -1
        self._popup: AcrylicPanel | None = None
        self._before_popup = None
        self._selection_changed = None
        self.clicked.connect(self._show_popup)

    def addItem(self, text: str, user_data=None) -> None:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        self._items.append((text, user_data))
        if self._current_index < 0:
            self.setCurrentIndex(0)

    def clear(self) -> None:
        self._items.clear()
        self._current_index = -1
        self.setText(self._placeholder)
        self._close_popup()

    def count(self) -> int:
        return len(self._items)

    def currentText(self) -> str:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        if 0 <= self._current_index < len(self._items):
            return self._items[self._current_index][0]
        return ""

    def currentData(self):  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        if 0 <= self._current_index < len(self._items):
            return self._items[self._current_index][1]
        return None

    def findText(self, text: str) -> int:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        for index, (item_text, _data) in enumerate(self._items):
            if item_text == text:
                return index
        return -1

    def findData(self, data) -> int:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        for index, (_text, item_data) in enumerate(self._items):
            if item_data == data:
                return index
        return -1

    def setCurrentIndex(self, index: int) -> None:  # noqa: N802 - 保持类似 QComboBox 的调用方式。
        if 0 <= index < len(self._items):
            self._current_index = index
            self.setText(self._items[index][0])
        else:
            self._current_index = -1
            self.setText(self._placeholder)

    def set_items(self, items: list[str]) -> None:
        self.clear()
        for item in items:
            self.addItem(item, item)

    def set_before_popup_callback(self, callback) -> None:
        self._before_popup = callback

    def set_selection_changed_callback(self, callback) -> None:
        self._selection_changed = callback

    def setEnabled(self, enabled: bool) -> None:  # noqa: N802 - Qt 接口沿用 Qt 命名。
        if not enabled:
            self._close_popup()
        super().setEnabled(enabled)

    def _show_popup(self) -> None:
        if not self.isEnabled():
            return
        if self._popup is not None:
            self._close_popup()
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
            for index, (text, _data) in enumerate(self._items):
                option = AcrylicTextButton(text)
                option.setObjectName("selectOption")
                option.setFixedHeight(30)
                option.clicked.connect(lambda _checked=False, item_index=index: self._select_item(item_index))
                layout.addWidget(option)
        else:
            empty_label = QLabel("暂无可选项")
            empty_label.setObjectName("selectEmpty")
            empty_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(empty_label)

        popup.setFixedWidth(self.width())
        popup.adjustSize()
        popup.move(root.mapFromGlobal(self.mapToGlobal(self.rect().bottomLeft())))
        popup.raise_()
        self._popup = popup
        popup.show()

    def _select_item(self, index: int) -> None:
        self.setCurrentIndex(index)
        if self._selection_changed is not None:
            self._selection_changed(self.currentText())
        self._close_popup()

    def _close_popup(self) -> None:
        if self._popup is not None:
            self._popup.deleteLater()
            self._popup = None


class ChannelCard(AnimatedButtonMixin, QPushButton):
    """分光模块单个通道显示卡片，绘制风格沿用参考 LedCard。"""

    def __init__(self, channel_index: int) -> None:
        super().__init__()
        self.channel_index = channel_index
        self.setObjectName("ledCard")
        self.setCursor(Qt.ArrowCursor)
        self.setCheckable(False)
        self.setMinimumSize(132, 104)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._hovered = False
        self._state = "idle"
        self._value: float | None = None
        self._baseline: float | None = None
        self._stable: float | None = None
        self._shadow = apply_shadow(self, blur_radius=18, y_offset=4, alpha=22)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self._init_button_animation()
        self._refresh_shadow()

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
        blur_radius = _lerp(22, 27, self._hover_progress)
        y_offset = _lerp(5, 7, self._hover_progress)
        alpha = _lerp(24, 38, self._hover_progress)
        self._shadow.setBlurRadius(blur_radius)
        self._shadow.setOffset(0, y_offset)
        self._shadow.setColor(QColor(15, 23, 42, int(alpha)))
        self.update()

    def set_state(self, state: str) -> None:
        self._state = state
        self.update()

    def set_value(self, value: float | None) -> None:
        self._value = value
        self.update()

    def set_baseline(self, value: float | None) -> None:
        self._baseline = value
        self.update()

    def set_stable(self, value: float | None) -> None:
        self._stable = value
        self.update()

    def reset_measurement(self) -> None:
        self._state = "idle"
        self._value = None
        self._baseline = None
        self._stable = None
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt 事件函数沿用 Qt 命名。
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        text_rect = QRectF(self.rect()).adjusted(1.0, 1.0, -1.0, -1.0)
        rect = self._animated_rect(text_rect, scale_ratio=0.005)
        path = QPainterPath()
        path.addRoundedRect(rect, 18, 18)

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
        painter.setPen(QPen(border_color, 1.0))
        painter.drawPath(path)

        dot_color = DOT_COLORS.get(self._state, DOT_COLORS["idle"])
        painter.setPen(QPen(QColor(255, 255, 255, 210), 1.5))
        painter.setBrush(dot_color)
        painter.drawEllipse(QRectF(text_rect.right() - 24, text_rect.top() + 12, 9, 9))

        title_font = QFont(self.font())
        title_font.setPointSize(13)
        title_font.setWeight(QFont.Weight.DemiBold)
        value_font = QFont(self.font())
        value_font.setPointSize(14)
        value_font.setWeight(QFont.Weight.Normal)
        info_font = QFont(self.font())
        info_font.setPointSize(9)
        info_font.setWeight(QFont.Weight.Normal)

        value_text = "--" if self._value is None else f"{self._value:.6f}"
        baseline_text = "--" if self._baseline is None else f"{self._baseline:.6f}"
        stable_text = "--" if self._stable is None else f"{self._stable:.6f}"

        painter.setPen(QColor("#0f172a"))
        painter.setFont(title_font)
        painter.drawText(text_rect.adjusted(12, 12, -32, -text_rect.height() * 0.62), Qt.AlignLeft | Qt.AlignTop, f"CH{self.channel_index}")

        painter.setFont(value_font)
        painter.drawText(text_rect.adjusted(12, text_rect.height() * 0.28, -12, -text_rect.height() * 0.42), Qt.AlignLeft | Qt.AlignTop, value_text)

        painter.setFont(info_font)
        painter.setPen(QColor("#475569"))
        painter.drawText(text_rect.adjusted(12, text_rect.height() * 0.58, -12, -24), Qt.AlignLeft | Qt.AlignTop, f"基底 {baseline_text}")
        painter.drawText(text_rect.adjusted(12, text_rect.height() * 0.76, -12, -8), Qt.AlignLeft | Qt.AlignTop, f"稳定 {stable_text}")
