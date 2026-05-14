from __future__ import annotations

import time

from PySide6.QtCore import QThread, Signal

from led_debugger.protocol import (
    ADC_PAYLOAD_LENGTH,
    RESPONSE_HEADER,
    ProtocolError,
    build_read_adc_request,
    parse_adc_response,
)
from led_debugger.settings import (
    DEFAULT_BAUD_RATE,
    DEFAULT_READ_TIMEOUT_SECONDS,
    DEFAULT_ROUND_DELAY_SECONDS,
)

import serial
from serial.tools import list_ports


BAUD_RATE = DEFAULT_BAUD_RATE
POLL_CHANNELS = range(1, 17)
READ_TIMEOUT_SECONDS = DEFAULT_READ_TIMEOUT_SECONDS
ROUND_DELAY_SECONDS = DEFAULT_ROUND_DELAY_SECONDS


def list_serial_ports() -> list[str]:
    """返回当前 Windows 可见串口名称，例如 COM3。"""
    return [port.device for port in list_ports.comports()]


class SerialPoller(QThread):
    """后台串口轮询线程。

    线程连接硬件后按 1-16 通道循环发送读取命令，并通过信号把每个通道的
    最新 16 路电流值推送给界面，避免串口等待阻塞 GUI。
    """

    reading_received = Signal(int, object)
    status_changed = Signal(str)
    error_occurred = Signal(str)

    def __init__(
        self,
        port_name: str,
        parent=None,
        read_timeout_seconds: float = READ_TIMEOUT_SECONDS,
        round_delay_seconds: float = ROUND_DELAY_SECONDS,
    ) -> None:
        super().__init__(parent)
        self.port_name = port_name
        self.read_timeout_seconds = read_timeout_seconds
        self.round_delay_seconds = round_delay_seconds
        self._running = False

    def stop(self) -> None:
        """请求线程停止；实际关闭串口在 run 内完成。"""
        self._running = False

    def run(self) -> None:
        self._running = True
        try:
            with serial.Serial(
                self.port_name,
                baudrate=BAUD_RATE,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE,
                timeout=self.read_timeout_seconds,
                write_timeout=self.read_timeout_seconds,
            ) as port:
                self.status_changed.emit(f"状态：已连接 {self.port_name}")
                self._poll_loop(port)
        except Exception as exc:
            if self._running:
                self.error_occurred.emit(f"串口错误：{exc}")
        finally:
            self._running = False
            self.status_changed.emit("状态：未连接")

    def _poll_loop(self, port) -> None:
        """持续轮询 16 个通道，切换页面时界面可直接显示最新缓存值。"""
        while self._running:
            for channel in POLL_CHANNELS:
                if not self._running:
                    break
                try:
                    port.reset_input_buffer()
                    port.write(build_read_adc_request(channel))
                    frame = self._read_response_frame(port)
                    reading = parse_adc_response(frame, expected_channel=channel)
                    self.reading_received.emit(reading.channel, reading.currents_ma)
                except (TimeoutError, ProtocolError) as exc:
                    self.error_occurred.emit(f"通道 {channel} 读取失败：{exc}")
                except Exception as exc:
                    self.error_occurred.emit(f"通道 {channel} 串口异常：{exc}")
                    return
            time.sleep(self.round_delay_seconds)

    def _read_response_frame(self, port) -> bytes:
        """同步读取一帧 5A A5 开头的 ADC 返回数据。"""
        header = self._read_until_header(port)
        rest_length = 1 + 1 + 1 + ADC_PAYLOAD_LENGTH + 1
        rest = port.read(rest_length)
        if len(rest) != rest_length:
            raise TimeoutError(f"返回帧长度不足：收到 {len(rest)} / {rest_length} 字节")
        return header + rest

    def _read_until_header(self, port) -> bytes:
        """从串口流中查找返回帧头，丢弃前面的杂散字节。"""
        deadline = time.monotonic() + self.read_timeout_seconds
        buffer = bytearray()
        while time.monotonic() < deadline and self._running:
            chunk = port.read(1)
            if not chunk:
                continue
            buffer.extend(chunk)
            if len(buffer) > len(RESPONSE_HEADER):
                del buffer[0 : len(buffer) - len(RESPONSE_HEADER)]
            if bytes(buffer) == RESPONSE_HEADER:
                return RESPONSE_HEADER
        raise TimeoutError("未收到返回帧头 5A A5")
