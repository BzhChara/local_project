from __future__ import annotations

import struct
from dataclasses import dataclass


REQUEST_HEADER = bytes((0xA5, 0x5A))
RESPONSE_HEADER = bytes((0x5A, 0xA5))
FUNCTION_LED_CURRENT_ADC = 0x02
CHANNEL_MIN = 0x01
CHANNEL_MAX = 0x10
ADC_VALUE_COUNT = 16
ADC_PAYLOAD_LENGTH = ADC_VALUE_COUNT * 4
# 协议字段命名为 IHH/IHL/ILH/ILL，按“高字节在前”默认解析为大端 float。
FLOAT_ENDIAN = ">"


class ProtocolError(ValueError):
    """Raised when a device frame does not match the expected protocol."""


@dataclass(frozen=True)
class ChannelReading:
    """单个通道的一次采集结果。"""

    channel: int
    currents_ma: tuple[float, ...]


def checksum(frame_without_checksum: bytes) -> int:
    """SUM 校验：从帧头开始逐字节累加，取低 8 位。"""
    return sum(frame_without_checksum) & 0xFF


def build_read_adc_request(channel: int) -> bytes:
    """生成读取指定通道 LED 电流 ADC 的请求帧。

    帧格式来自协议表：
    A5 5A XX 02 02 00 00 SUM
    """
    _validate_channel(channel)
    frame = bytes(
        (
            *REQUEST_HEADER,
            channel,
            FUNCTION_LED_CURRENT_ADC,
            0x02,
            0x00,
            0x00,
        )
    )
    return frame + bytes((checksum(frame),))


def parse_adc_response(frame: bytes | bytearray, expected_channel: int | None = None) -> ChannelReading:
    """解析主板返回的 LED 电流 ADC 数据帧。

    返回帧格式：
    5A A5 XX 02 40 + 16 个大端 float + SUM
    """
    raw = bytes(frame)
    if expected_channel is not None:
        _validate_channel(expected_channel)

    # 最小长度包含帧头、通道、功能码、长度和校验位。
    minimum_length = 2 + 1 + 1 + 1 + 1
    if len(raw) < minimum_length:
        raise ProtocolError(f"Frame too short: {len(raw)} bytes")

    if raw[:2] != RESPONSE_HEADER:
        raise ProtocolError("Invalid response header")

    # 通道范围固定为 0x01 到 0x10，和协议表 XX：0x01....0x10 保持一致。
    channel = raw[2]
    _validate_channel(channel)
    if expected_channel is not None and channel != expected_channel:
        raise ProtocolError(f"Unexpected channel: got {channel}, expected {expected_channel}")

    function_code = raw[3]
    if function_code != FUNCTION_LED_CURRENT_ADC:
        raise ProtocolError(f"Unexpected function code: 0x{function_code:02X}")

    payload_length = raw[4]
    if payload_length != ADC_PAYLOAD_LENGTH:
        raise ProtocolError(f"Unexpected payload length: {payload_length}")

    expected_length = 2 + 1 + 1 + 1 + payload_length + 1
    if len(raw) != expected_length:
        raise ProtocolError(f"Unexpected frame length: got {len(raw)}, expected {expected_length}")

    # 先校验 SUM，再解析 payload，避免错误帧进入后续计算链路。
    expected_checksum = checksum(raw[:-1])
    actual_checksum = raw[-1]
    if actual_checksum != expected_checksum:
        raise ProtocolError(
            f"Invalid checksum: got 0x{actual_checksum:02X}, expected 0x{expected_checksum:02X}"
        )

    payload_start = 5
    payload_end = payload_start + payload_length
    # 每个 LED 电流值占 4 字节，当前按 IEEE754 大端 float 解析，单位暂定 mA。
    values = struct.unpack(f"{FLOAT_ENDIAN}{ADC_VALUE_COUNT}f", raw[payload_start:payload_end])
    return ChannelReading(channel=channel, currents_ma=tuple(values))


def _validate_channel(channel: int) -> None:
    """校验协议通道号范围。"""
    if not CHANNEL_MIN <= channel <= CHANNEL_MAX:
        raise ProtocolError(f"Channel out of range: {channel}")
