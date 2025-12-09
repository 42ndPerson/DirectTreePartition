cdef class Point:
    cdef public double x
    cdef public double y
    cpdef tuple tuple(self)
    
    @staticmethod
    cdef double dist(Point p1, Point p2)
    
    @staticmethod
    cdef double src_dest_rad(Point src, Point dest)
    
    @staticmethod
    cdef double normalized_angle_diff(double angle_from, double angle_to)
