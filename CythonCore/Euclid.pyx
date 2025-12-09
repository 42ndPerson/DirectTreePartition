# cython: language_level=3
from libc.math cimport sqrt, atan2, pi

cdef class Point:
    def __init__(self, double x, double y):
        self.x = x
        self.y = y

    cpdef tuple tuple(self):
        return (self.x, self.y)

    @staticmethod
    cdef double dist(Point p1, Point p2):
        return sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)
    
    @staticmethod
    cdef double src_dest_rad(Point src, Point dest):
        return atan2(dest.y - src.y, dest.x - src.x)
    
    @staticmethod
    cdef double normalized_angle_diff(double angle_from, double angle_to):
        cdef double delta = angle_to - angle_from
        delta = (delta + pi) % (2 * pi) - pi
        return delta
