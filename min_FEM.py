import numpy as np

# minimal implementation with One square (four nodes)
"""
N4 —— N3
|     |
|     |
|     |
N1 —— N2

First triangle: globals N1, N2, N4
Second triangle: globals N2, N3, N4 -- n1, n2, n3 local

"""


class Node:
    def __init__(self, id, x, y):
        self.id = id  # global
        self.x = x
        self.y = y


class Triangle:  # element class
    def __init__(self, id, n1, n2, n3, k, h):
        self.id = id
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.k = k
        self.h = h
        self.centroidX = (n1.x + n2.x + n3.x) / 3
        self.centroidY = (n1.y + n2.y + n3.y) / 3
        self.gradient = None
        self.flux = None

    def node_ids(self):
        return np.array(
            [self.n1.id - 1, self.n2.id - 1, self.n3.id - 1]
        )  # 0-indexed for np arrays

    def area(self):
        x1, y1 = self.n1.x, self.n1.y
        x2, y2 = self.n2.x, self.n2.y
        x3, y3 = self.n3.x, self.n3.y
        return 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

    def a(self, i, j):
        return i.x * j.y - j.x * i.y

    def b(self, i, j):
        """X coord diff between each node of the triangle"""
        return i.y - j.y

    def c(self, i, j):
        """Y coord diff between each node of the triangle"""
        return i.x - j.x

    def a_coeffs(self):
        """Area coefficients for the triangle"""
        return [
            self.a(self.n2, self.n3),
            self.a(self.n3, self.n1),
            self.a(self.n1, self.n2),
        ]

    def b_coeffs(self):
        """b coefficients for the triangle"""
        return [
            self.b(self.n2, self.n3),
            self.b(self.n3, self.n1),
            self.b(self.n1, self.n2),
        ]

    def c_coeffs(self):
        """c coefficients for the triangle"""
        return [
            self.c(self.n3, self.n2),
            self.c(self.n1, self.n3),
            self.c(self.n2, self.n1),
        ]
