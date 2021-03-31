from math import ceil, log2
import numpy as np

"""
van Emde Boas Tree

Provides O(log(log(u)) query time for operations like INSERT, SEARCH, DELETE, SUCCESSOR, and PREDECESSOR
"""

class VEB:

    def high(self, x):
        return x >> (self.w // 2)

    def low(self, x):
        return x & (1 << (self.w // 2)) - 1

    def index(self, i, j):
        return i << (self.w // 2) | j

    def __init__(self, u):
        self.w = ceil(log2(u))
        self.min = self.max = None
        self.start = self.stop = None

        if self.w >= 1:
            self.cluster = {}
            self.summary = None

    # def member(self, x):
    #
    #     if x == self.min or x == self.max:
    #         return True
    #     elif self.w == 1:
    #         return False
    #     else:
    #         c = self.high(x)
    #         i = self.low(x)
    #         if c in self.cluster:
    #             return self.cluster[c].member(i)
    #         else:
    #             return False
    def member(self, x):

        if x == self.min or x == self.max:
            return (self.start, self.stop)
        elif self.w == 1:
            return False
        else:
            c = self.high(x)
            i = self.low(x)
            if c in self.cluster:
                return self.cluster[c].member(i)
            else:
                return False

    def insert(self, x, start, stop):
        if self.min is None:
            self.min = x
            self.max = x
            self.start = start
            self.stop = stop
            return
        else:
            if x < self.min:
                x, self.min = self.min, x
            c = self.high(x)
            i = self.low(x)
            if self.w > 1:
                if c not in self.cluster:
                    self.cluster[c] = VEB(2 ** (self.w // 2))
                if self.cluster[c].min is None:
                    if self.summary is None:
                        self.summary = VEB(2 ** (self.w // 2))
                    self.summary.insert(c, start, stop)
                if c not in self.cluster:
                    self.cluster[c] = VEB(2 ** (self.w // 2))

                self.cluster[c].insert(i, start, stop)

            if x > self.max:
                self.max = x
