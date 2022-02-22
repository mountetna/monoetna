import time


class RetryGater:
    def __init__(self, max: int, f=1, exp=2):
        self.max = max
        self.i = 0
        self.last = 0
        self.exp = exp
        self.f = f

    def reset(self):
        self.i = 0

    def gate(self, exception: Exception):
        self.i += 1
        if self.i > self.max:
            raise exception

        time.sleep((self.exp ** (self.i - 1)) * self.f)
