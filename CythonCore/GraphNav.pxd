from CythonCore.TwinGraph cimport TwinGraph, Vert, QuadEdge, EdgeDir, VertRole, reverse_dir
from CythonCore.RegionTree cimport RegionTree, Region, Edge

cdef class GraphNav

cpdef enum StartSelectionMethod:
    UNIFORM = 1
    CENTRAL = 2
    LINING = 3
    WALK = 4

cpdef enum MultiWalkStartBehavior:
    RESTART = 1
    HIT_POINT = 2

cdef class GraphNav:
    cdef public TwinGraph graph
    cdef public RegionTree region_tree
    cdef public object target_regions # deque[Region]
    cdef public set tree_verts
    cdef public set tree_edges
    cdef public set bridge_edges
    cdef public StartSelectionMethod start_selection_method
    cdef public MultiWalkStartBehavior multi_walk_start_behavior
    cdef public int multi_walk_attempts
    cdef public bint unwind_loops
    cdef public bint animating
    cdef public list animation_tracks
    
    cpdef list run_two_split_attempt(self)
    cpdef set walk_division_from(self, Region region_in, Vert origin_vert)
    cpdef set get_enclosing_perimeter_verts(self, Vert vert)
    cpdef set get_enclosed_primal_verts(self, list perim)
    cpdef tuple loop_erased_random_walk_from(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region)
    cpdef tuple _loop_erased_random_walk_from_lazy(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region)
    cpdef tuple _loop_erased_random_walk_from_unwind(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region)
    cpdef QuadEdge commit_walk(self, list walk_edges, set consumed_verts)
    cpdef list develop_region(self, Region src_region, QuadEdge start_edge, EdgeDir start_dir)
    cpdef list traverse_clockwise_loop(self, QuadEdge start_edge, EdgeDir edge_dir, set bridge_boundaries)
