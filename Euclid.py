import math

class Euclid:
    class Point:
        x: float
        y: float

        def __init__(self, x: float, y: float) -> None:
            self.x = x
            self.y = y

        @staticmethod
        def dist(p1: Euclid.Point, p2: Euclid.Point) -> float:
            return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)
        
        @staticmethod
        def srcDestRad(src: Euclid.Point, dest: Euclid.Point) -> float:
            return math.atan2(dest.y - src.y, dest.x - src.x)