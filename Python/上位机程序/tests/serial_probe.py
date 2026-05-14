import struct
import time

import serial

PORT = "COM5"
BAUD = 115200
TIMEOUT = 0.4

REQ_HEADER = bytes((0xA5, 0x5A))
RESP_HEADER = bytes((0x5A, 0xA5))
FUNC = 0x02
PAYLOAD_LEN = 64


def checksum(data: bytes) -> int:
    return sum(data) & 0xFF


def build_req(ch: int) -> bytes:
    frame = bytes((*REQ_HEADER, ch, FUNC, 0x02, 0x00, 0x00))
    return frame + bytes((checksum(frame),))


def read_until_header(port):
    deadline = time.monotonic() + TIMEOUT
    buf = bytearray()
    while time.monotonic() < deadline:
        b = port.read(1)
        if not b:
            continue
        buf.extend(b)
        if len(buf) > 2:
            del buf[0 : len(buf) - 2]
        if bytes(buf) == RESP_HEADER:
            return bytes(buf)
    raise TimeoutError("no 5A A5 header")


def parse_frame(frame: bytes):
    if len(frame) != 70:
        return None, f"bad length {len(frame)}"
    if frame[:2] != RESP_HEADER:
        return None, "bad header"

    expected = checksum(frame[:-1])
    actual = frame[-1]
    if expected != actual:
        return None, f"bad checksum actual=0x{actual:02X} expected=0x{expected:02X}"

    if frame[3] != FUNC or frame[4] != PAYLOAD_LEN:
        return None, f"bad func/len func=0x{frame[3]:02X} len={frame[4]}"

    vals = struct.unpack(">16f", frame[5:69])
    return vals, "ok"


with serial.Serial(PORT, BAUD, timeout=TIMEOUT, write_timeout=TIMEOUT) as port:
    print(f"OPEN {PORT} {BAUD}")
    for ch in range(1, 17):
        try:
            port.reset_input_buffer()
            req = build_req(ch)
            port.write(req)

            header = read_until_header(port)
            rest = port.read(68)
            frame = header + rest

            vals, status = parse_frame(frame)
            print(f"CH{ch:02d} REQ={req.hex(' ').upper()} LEN={len(frame)} STATUS={status}")
            print(f"CH{ch:02d} RAW={frame.hex(' ').upper()}")

            if vals is not None:
                formatted = ", ".join(f"{v:.6g}" for v in vals)
                nonzero = [i + 1 for i, v in enumerate(vals) if abs(v) >= 1e-9]
                print(f"CH{ch:02d} VALUES=[{formatted}] NONZERO={nonzero}")

        except Exception as exc:
            print(f"CH{ch:02d} ERROR={type(exc).__name__}: {exc}")

        time.sleep(0.02)