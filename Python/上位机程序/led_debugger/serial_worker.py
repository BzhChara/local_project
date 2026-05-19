from __future__ import annotations

import time
from dataclasses import dataclass

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
    DEFAULT_CHANNEL_DELAY_SECONDS,
    DEFAULT_READ_TIMEOUT_SECONDS,
    DEFAULT_ROUND_DELAY_SECONDS,
)

import serial
from serial.tools import list_ports


BAUD_RATE = DEFAULT_BAUD_RATE
POLL_CHANNELS = range(1, 17)
READ_TIMEOUT_SECONDS = DEFAULT_READ_TIMEOUT_SECONDS
ROUND_DELAY_SECONDS = DEFAULT_ROUND_DELAY_SECONDS
CHANNEL_DELAY_SECONDS = DEFAULT_CHANNEL_DELAY_SECONDS


@dataclass(frozen=True)
class RawSerialRecord:
    """单个通道最近一次串口原始返回。"""

    channel: int
    timestamp: float
    status: str
    raw_bytes: bytes


class SerialReadTimeout(TimeoutError):
    """读取串口返回帧失败，并携带实际收到的字节。"""

    def __init__(self, message: str, raw_bytes: bytes) -> None:
        super().__init__(message)
        self.raw_bytes = raw_bytes


def list_serial_ports() -> list[str]:
    """返回当前 Windows 可见串口名称，例如 COM3。"""
    return [port.device for port in list_ports.comports()]


class SerialPoller(QThread):
    """后台串口轮询线程。

    线程连接硬件后按 1-16 通道循环发送读取命令，并通过信号把每个通道的
    最新 16 路电流值推送给界面，避免串口等待阻塞 GUI。
    """

    reading_received = Signal(int, object)
    raw_response_received = Signal(int, object)
    status_changed = Signal(str)
    error_occurred = Signal(str)

    def __init__(
        self,
        port_name: str,
        parent=None,
        read_timeout_seconds: float = READ_TIMEOUT_SECONDS,
        round_delay_seconds: float = ROUND_DELAY_SECONDS,
        channel_delay_seconds: float = CHANNEL_DELAY_SECONDS,
    ) -> None:
        super().__init__(parent)
        self.port_name = port_name
        self.read_timeout_seconds = read_timeout_seconds
        self.round_delay_seconds = round_delay_seconds
        self.channel_delay_seconds = channel_delay_seconds
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
                    frame, raw_bytes = self._read_response_frame(port)
                    try:
                        reading = parse_adc_response(frame, expected_channel=channel)
                    except ProtocolError as exc:
                        self._emit_raw_response(channel, f"解析失败：{exc}", raw_bytes)
                        raise
                    self._emit_raw_response(channel, "解析成功", raw_bytes)
                    self.reading_received.emit(reading.channel, reading.currents_ma)
                except TimeoutError as exc:
                    self._emit_raw_response(channel, str(exc), getattr(exc, "raw_bytes", b""))
                    self.error_occurred.emit(f"通道 {channel} 读取失败：{exc}")
                except ProtocolError as exc:
                    self.error_occurred.emit(f"通道 {channel} 读取失败：{exc}")
                except Exception as exc:
                    self.error_occurred.emit(f"通道 {channel} 串口异常：{exc}")
                    return
                if self._running and self.channel_delay_seconds > 0:
                    time.sleep(self.channel_delay_seconds)
            time.sleep(self.round_delay_seconds)

    def _read_response_frame(self, port) -> tuple[bytes, bytes]:
        """同步读取一帧 5A A5 开头的 ADC 返回数据。"""
        header, captured = self._read_until_header(port)
        rest_length = 1 + 1 + 1 + ADC_PAYLOAD_LENGTH + 1
        rest = port.read(rest_length)
        raw_bytes = captured + rest
        if len(rest) != rest_length:
            raise SerialReadTimeout(f"返回帧长度不足：收到 {len(rest)} / {rest_length} 字节", raw_bytes)
        return header + rest, raw_bytes

    def _read_until_header(self, port) -> tuple[bytes, bytes]:
        """从串口流中查找返回帧头，丢弃前面的杂散字节。"""
        deadline = time.monotonic() + self.read_timeout_seconds
        buffer = bytearray()
        captured = bytearray()
        while time.monotonic() < deadline and self._running:
            chunk = port.read(1)
            if not chunk:
                continue
            captured.extend(chunk)
            buffer.extend(chunk)
            if len(buffer) > len(RESPONSE_HEADER):
                del buffer[0 : len(buffer) - len(RESPONSE_HEADER)]
            if bytes(buffer) == RESPONSE_HEADER:
                return RESPONSE_HEADER, bytes(captured)
        raise SerialReadTimeout("未收到返回帧头 5A A5", bytes(captured))

    def _emit_raw_response(self, channel: int, status: str, raw_bytes: bytes) -> None:
        """把诊断用原始返回发给界面。"""
        self.raw_response_received.emit(
            channel,
            RawSerialRecord(
                channel=channel,
                timestamp=time.time(),
                status=status,
                raw_bytes=bytes(raw_bytes),
            ),
        )
