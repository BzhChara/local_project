from __future__ import annotations

from datetime import datetime

from PySide6.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout

from led_debugger.serial_worker import RawSerialRecord
from led_debugger.widgets import AcrylicPanel


class RawSerialPanel(AcrylicPanel):
    """右侧串口原始返回诊断面板。"""

    def __init__(self) -> None:
        super().__init__(radius=20, top_alpha=118, bottom_alpha=88)
        self.current_channel = 1
        self.current_record: RawSerialRecord | None = None
        self.setObjectName("rawSerialPanel")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        self.title = QLabel("通道 1 原始返回")
        self.title.setObjectName("curveTitle")

        self.text_box = QPlainTextEdit()
        self.text_box.setObjectName("rawSerialText")
        self.text_box.setReadOnly(True)
        self.text_box.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.text_box.setPlainText("暂无该通道返回数据")

        layout.addWidget(self.title)
        layout.addWidget(self.text_box, 1)

    def set_channel(self, channel: int) -> None:
        """切换当前显示通道。"""
        self.current_channel = channel
        self.title.setText(f"通道 {channel} 原始返回")
        self._render()

    def set_record(self, record: RawSerialRecord | None) -> None:
        """更新当前通道最近一次原始返回。"""
        self.current_record = record
        self._render()

    def _render(self) -> None:
        if self.current_record is None:
            self.text_box.setPlainText("暂无该通道返回数据")
            return

        record = self.current_record
        timestamp_text = datetime.fromtimestamp(record.timestamp).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
        hex_text = format_hex(record.raw_bytes)
        self.text_box.setPlainText(
            "\n".join(
                (
                    f"时间：{timestamp_text}",
                    f"状态：{record.status}",
                    f"长度：{len(record.raw_bytes)} 字节",
                    "RX：",
                    hex_text if hex_text else "无返回",
                )
            )
        )


def format_hex(raw_bytes: bytes, bytes_per_line: int = 16) -> str:
    """把原始字节格式化为多行 HEX。"""
    if not raw_bytes:
        return ""
    groups = []
    for start in range(0, len(raw_bytes), bytes_per_line):
        chunk = raw_bytes[start : start + bytes_per_line]
        groups.append(" ".join(f"{byte:02X}" for byte in chunk))
    return "\n".join(groups)
