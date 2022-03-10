import io
from typing import Iterable


class IterStream(io.RawIOBase):
    def __init__(self, iterable):
        self.iterable = iterable
        self.leftover = None

    def readable(self):
        return True

    def readinto(self, b):
        try:
            l = len(b)  # We're supposed to return at most this much
            chunk = self.leftover or next(self.iterable)
            output, self.leftover = chunk[:l], chunk[l:]
            b[: len(output)] = output
            return len(output)
        except StopIteration:
            return 0  # indicate EOF


def iterable_to_stream(
    iterable: Iterable[bytes], buffer_size=io.DEFAULT_BUFFER_SIZE
) -> io.BufferedReader:
    """
    Lets you use an iterable (e.g. a generator) that yields bytestrings as a read-only
    input stream.

    The stream implements Python 3's newer I/O API (available in Python 2's io module).
    For efficiency, the stream is buffered.
    """
    return io.BufferedReader(IterStream(iterable), buffer_size=buffer_size)
