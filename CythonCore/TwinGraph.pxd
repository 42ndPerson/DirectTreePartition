from CythonCore.Euclid cimport Point

# Forward declarations
cdef class Vert
cdef class QuadEdge
cdef class MapTree
cdef class MapTreeVert
cdef class MapTreeEdge
cdef class TwinGraph

cpdef enum EdgeDir:
    AB = 1
    BA = 2

cpdef EdgeDir reverse_dir(EdgeDir dir)

cpdef enum VertRole:
    PRIMAL = 1
    DUAL = -1
    DUAL_EXTERIOR = -2

cdef class Vert:
    cdef public Point point
    cdef public int weight
    cdef public VertRole role
    cdef public list cc_edges # List[QuadEdge]
    cdef public MapTreeVert map_tree_vert
    cdef public int index
    cdef public str id_str
    
    cpdef void register_edges(self, list edges)
    cpdef void link_edges(self)

cdef class QuadEdge:
    cdef public Vert primal_A
    cdef public Vert primal_B
    cdef public QuadEdge primal_A_cc_next
    cdef public QuadEdge primal_A_cc_prev
    cdef public QuadEdge primal_B_cc_next
    cdef public QuadEdge primal_B_cc_prev
    
    cdef public Vert dual_AB
    cdef public Vert dual_BA
    cdef public QuadEdge dual_AB_cc_next
    cdef public QuadEdge dual_AB_cc_prev
    cdef public QuadEdge dual_BA_cc_next
    cdef public QuadEdge dual_BA_cc_prev
    
    cdef public object dual_AB_annotation # Optional[int]
    cdef public MapTreeEdge map_tree_edge
    cdef public str id_str
    
    cpdef tuple get_primal_vert_pair(self, EdgeDir dir)
    cpdef tuple get_dual_vert_pair(self, EdgeDir dir)
    cpdef tuple get_primal_dest_from(self, Vert src)
    cpdef tuple get_dual_dest_from(self, Vert src)
    cpdef tuple get_primal_rad_from(self, Vert src)
    cpdef tuple get_dual_rad_from(self, Vert src)
    cpdef double get_primal_rad_along(self, EdgeDir dir)
    cpdef double get_dual_rad_along(self, EdgeDir dir)
    cpdef tuple get_primal_cc_next_edge(self, Vert vert)
    cpdef tuple get_dual_cc_next_edge(self, Vert vert)
    cpdef tuple get_primal_cc_prev_edge(self, Vert vert)
    cpdef tuple get_dual_cc_prev_edge(self, Vert vert)
    cpdef Vert get_common_primal_vert(self, QuadEdge other)
    cpdef Vert get_common_dual_vert(self, QuadEdge other)

cdef class MapTreeVert:
    cdef public Vert twin_vert
    cdef public int dfs_in
    cdef public int dfs_out
    cdef public list out_edges # List[MapTreeEdge]
    cdef public MapTreeEdge in_edge
    
    cpdef MapTreeEdge add_child(self, MapTreeVert child, QuadEdge twin_edge, EdgeDir edge_dir)

cdef class MapTreeEdge:
    cdef public MapTreeVert parent
    cdef public MapTreeVert child
    cdef public QuadEdge twin_edge
    cdef public EdgeDir edge_dir
    
    cpdef MapTreeVert get_dest_from(self, MapTreeVert src)

cdef class MapTree:
    cdef public MapTreeVert root
    cdef public set verts
    cdef public set edges
    cdef public set quad_edges
    
    cpdef void collect_tree(self, MapTreeVert vert)

cdef class TwinGraph:
    cdef public set primalVerts
    cdef public int total_primal_weight
    cdef public set dualVerts
    cdef public Vert external_dual_vert
    cdef public set edges
    cdef public double lowerXBound
    cdef public double upperXBound
    cdef public double lowerYBound
    cdef public double upperYBound
    cdef public MapTree primalMap
    cdef public MapTree dualMap
    cdef public list idx_to_dfs_in
    cdef public list idx_to_dfs_out
    cdef public bint animating
    cdef public list animation_tracks
    cdef public list dual_verts_list
    
    cpdef void construct_dual(self)
    cpdef void construct_face_loop(self, QuadEdge src_edge, EdgeDir src_dir, set edges_AB, set edges_BA)
    cpdef MapTree generate_map_tree(self, Vert root_vert, VertRole role)
    cdef int _dfs(self, Vert parent_vert, MapTreeVert parent_tree_vert, int index, set visited_verts, set visited_edges, set tree_verts, VertRole role)
    cpdef void populate_index_to_dfs_lookup(self)
    cpdef void annotate_dual_edges_for_primal_counting(self)
    cpdef int count_primal_verts_within_perim(self, list perim)
    cpdef list get_verts_within_radius(self, Point src, double radius, VertRole role)
    cpdef Vert get_closest_vert(self, Point src, VertRole role)
