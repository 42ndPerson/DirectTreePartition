from __future__ import annotations
from typing import Tuple

import math

class Point:
    x: float
    y: float

    def __init__(self, x: float, y: float) -> None:
        self.x = x
        self.y = y

    def tuple(self) -> Tuple[float, float]:
        return (self.x, self.y)

    @staticmethod
    def dist(p1: Point, p2: Point) -> float:
        return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)
    
    @staticmethod
    def src_dest_rad(src: Point, dest: Point) -> float:
        return math.atan2(dest.y - src.y, dest.x - src.x)