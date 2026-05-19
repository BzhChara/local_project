import pytest

from led_debugger.protocol import RESPONSE_HEADER
from led_debugger.serial_worker import SerialPoller, SerialReadTimeout


class FakePort:
    def __init__(self, chunks: list[bytes]) -> None:
        self.chunks = chunks

    def read(self, _length: int) -> bytes:
        if not self.chunks:
            return b""
        return self.chunks.pop(0)


def test_read_until_header_returns_full_captured_bytes() -> None:
    poller = SerialPoller("COM_TEST", read_timeout_seconds=0.01)
    poller._running = True

    header, captured = poller._read_until_header(FakePort([b"\x00", b"\x11", b"\x5A", b"\xA5"]))

    assert header == RESPONSE_HEADER
    assert captured == b"\x00\x11\x5A\xA5"


def test_read_until_header_timeout_keeps_received_bytes() -> None:
    poller = SerialPoller("COM_TEST", read_timeout_seconds=0.001)
    poller._running = True

    with pytest.raises(SerialReadTimeout) as exc_info:
        poller._read_until_header(FakePort([b"\x12", b"\x34"]))

    assert str(exc_info.value) == "未收到返回帧头 5A A5"
    assert exc_info.value.raw_bytes == b"\x12\x34"


def test_read_response_frame_timeout_keeps_partial_frame() -> None:
    poller = SerialPoller("COM_TEST", read_timeout_seconds=0.01)
    poller._running = True

    with pytest.raises(SerialReadTimeout) as exc_info:
        poller._read_response_frame(FakePort([b"\x5A", b"\xA5", b"\x01\x02"]))

    assert str(exc_info.value) == "返回帧长度不足：收到 2 / 68 字节"
    assert exc_info.value.raw_bytes == b"\x5A\xA5\x01\x02"
