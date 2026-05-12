import struct

import pytest

from led_debugger.protocol import (
    ADC_PAYLOAD_LENGTH,
    FUNCTION_LED_CURRENT_ADC,
    RESPONSE_HEADER,
    ProtocolError,
    build_read_adc_request,
    checksum,
    parse_adc_response,
)


def make_adc_response(channel: int, values: list[float]) -> bytes:
    """构造一帧合法的 ADC 返回帧，供协议解析测试使用。"""
    payload = struct.pack(">16f", *values)
    frame = bytes((*RESPONSE_HEADER, channel, FUNCTION_LED_CURRENT_ADC, ADC_PAYLOAD_LENGTH)) + payload
    return frame + bytes((checksum(frame),))


def test_checksum_uses_low_byte_of_sum() -> None:
    assert checksum(bytes.fromhex("A5 5A 01 02 02 00 00")) == 0x04


def test_build_read_adc_request() -> None:
    assert build_read_adc_request(1) == bytes.fromhex("A5 5A 01 02 02 00 00 04")
    assert build_read_adc_request(16) == bytes.fromhex("A5 5A 10 02 02 00 00 13")


def test_parse_adc_response_reads_big_endian_float_values() -> None:
    values = [float(i) for i in range(1, 17)]
    reading = parse_adc_response(make_adc_response(3, values), expected_channel=3)

    assert reading.channel == 3
    assert reading.currents_ma == pytest.approx(values)


def test_parse_adc_response_rejects_wrong_checksum() -> None:
    values = [float(i) for i in range(1, 17)]
    frame = bytearray(make_adc_response(1, values))
    frame[-1] ^= 0xFF

    with pytest.raises(ProtocolError, match="Invalid checksum"):
        parse_adc_response(frame)


def test_parse_adc_response_rejects_unexpected_channel() -> None:
    values = [float(i) for i in range(1, 17)]

    with pytest.raises(ProtocolError, match="Unexpected channel"):
        parse_adc_response(make_adc_response(2, values), expected_channel=1)


def test_channel_range_is_1_to_16() -> None:
    with pytest.raises(ProtocolError, match="Channel out of range"):
        build_read_adc_request(0)

    with pytest.raises(ProtocolError, match="Channel out of range"):
        build_read_adc_request(17)
