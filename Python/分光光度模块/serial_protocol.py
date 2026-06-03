import struct

import serial


SERIAL_BAUD = 115200
SERIAL_BYTESIZE = serial.EIGHTBITS
SERIAL_PARITY = serial.PARITY_NONE
SERIAL_STOPBITS = serial.STOPBITS_ONE


def checksum(frame_without_sum: bytes) -> int:
    return sum(frame_without_sum) & 0xFF


def build_cmd(channel: int, func: int, data: int) -> bytes:
    frame = bytearray([0xA5, 0x5A, channel & 0xFF, func & 0xFF, 0x02])
    frame.append((data >> 8) & 0xFF)
    frame.append(data & 0xFF)
    frame.append(checksum(frame))
    return bytes(frame)


def frame_to_hex(frame: bytes) -> str:
    return " ".join(f"{b:02X}" for b in frame)


def parse_frames(buffer: bytearray):
    frames = []
    while True:
        start = buffer.find(b"\x5A\xA5")
        if start < 0:
            if len(buffer) > 1:
                del buffer[:-1]
            return frames

        if start > 0:
            del buffer[:start]

        if len(buffer) < 6:
            return frames

        payload_len = buffer[4]
        total_len = 6 + payload_len
        if len(buffer) < total_len:
            return frames

        frame = bytes(buffer[:total_len])
        del buffer[:total_len]

        if checksum(frame[:-1]) != frame[-1]:
            print(f"丢弃校验错误帧: {frame_to_hex(frame)}")
            continue

        frames.append(frame)


def parse_absorbance_frame(frame: bytes):
    func = frame[3]
    payload_len = frame[4]
    if func not in (0x56, 0x57) or payload_len < 96:
        return None

    payload = frame[5:101]
    values = list(struct.unpack("<24f", payload))
    return [abs(v) for v in values]
