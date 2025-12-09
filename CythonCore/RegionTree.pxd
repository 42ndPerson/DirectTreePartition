from CythonCore.Euclid cimport Point
from CythonCore.TwinGraph cimport TwinGraph, Vert, QuadEdge, EdgeDir, VertRole, MapTreeVert

cdef class RegionTree
cdef class Region
cdef class Edge
cdef class BoundaryElement

cdef class BoundaryElement:
    cdef public int dfs_index
    cdef public int boundary_type
    cdef public Vert vert

cdef class RegionTree:
    cdef public set regions
    cdef public set edges
    cdef public TwinGraph graph
    cdef public double total_weight
    cdef public double epsilon
    cdef public set edge_centers
    cdef public str id_str
    
    cpdef void add_region(self, Region vert)
    cpdef void remove_region(self, Region vert)
    cpdef void add_edge(self, Edge edge)
    cpdef void remove_edge(self, Edge edge)
    cpdef void debug_inventory(self, Edge edge)

cdef class Region:
    cdef public int weight
    cdef public list dual_perimeter # List[Tuple[QuadEdge, EdgeDir]]
    cdef public set dual_perimeter_edges # Set[QuadEdge]
    cdef public set region_edges # Set[Edge]
    cdef public list cc_region_edges # List[Edge]
    cdef public set bridge_set # Set[QuadEdge]
    cdef public dict bridge_to_region_edge_map # Dict[QuadEdge, Edge]
    cdef public RegionTree region_tree
    cdef public Point point
    cdef public str id_str
    
    cpdef void register_edge(self, Edge edge)
    cpdef void unregister_edge(self, Edge edge)
    cpdef void generate_cc_region_edges(self)
    cpdef void calc_point(self)
    cpdef set get_perimeter_verts(self)
    cpdef tuple get_interior_lining_verts(self)
    cpdef Vert get_uniform_interior_vert(self)
    cpdef Vert get_central_region_vert(self)
    cpdef Vert get_walk_random_interior_vert(self, int steps=*)
    cdef Vert _find_interior_start_vert(self)
    cdef Vert _walk_interior_optimistic(self, Vert start, int step_cap, set perimeter_verts)
    cpdef bint check_center(self)
    cpdef bint validate_split_possible(self)

cdef class Edge:
    cdef public Region end_A
    cdef public Region end_B
    cdef public QuadEdge twin_graph_edge
    cdef public object ab_weight_differential # Optional[int]
    cdef public RegionTree region_tree
    cdef public str id_str
    
    cpdef tuple get_dest_from(self, Region src)
    cpdef object get_opposite_side_weight(self, Region src)
    cpdef void calculate_weight_differential_from_dir(self, EdgeDir dir)
