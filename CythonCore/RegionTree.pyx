# cython: language_level=3
import random
import bisect
from itertools import pairwise
from collections import deque
from libc.math cimport pi, sqrt, atan2

from CythonCore.Euclid cimport Point
from CythonCore.TwinGraph cimport TwinGraph, Vert, QuadEdge, EdgeDir, VertRole, MapTreeVert, reverse_dir
from operator import attrgetter

# Global counters
cdef int _region_counter = 0
cdef int _edge_counter = 0
cdef int _rtree_counter = 0

cdef class BoundaryElement:
    def __init__(self, int dfs_index, int boundary_type, Vert vert):
        self.dfs_index = dfs_index
        self.boundary_type = boundary_type
        self.vert = vert
    
    def __repr__(self):
        return f"BoundaryElement(dfs={self.dfs_index}, type={self.boundary_type}, vert={self.vert.id_str})"

cdef class RegionTree:
    def __init__(self, TwinGraph graph, double epsilon=0, list src_perimeter=None):
        global _rtree_counter
        self.graph = graph
        self.total_weight = graph.total_primal_weight
        self.epsilon = epsilon
        
        self.regions = set()
        self.add_region(Region(self.total_weight, src_perimeter if src_perimeter is not None else []))
        self.edges = set()
        self.edge_centers = set()
        
        self.id_str = f"RTree{_rtree_counter}"
        _rtree_counter += 1

    cpdef void add_region(self, Region vert):
        self.regions.add(vert)
        vert.region_tree = self

    cpdef void remove_region(self, Region vert):
        cdef Edge edge
        for edge in list(vert.region_edges):
            self.remove_edge(edge)
        self.regions.remove(vert)

    cpdef void add_edge(self, Edge edge):
        self.edges.add(edge)
        edge.end_A.register_edge(edge)
        edge.end_B.register_edge(edge)
        edge.region_tree = self

    cpdef void remove_edge(self, Edge edge):
        self.edges.remove(edge)
        edge.end_A.unregister_edge(edge)
        edge.end_B.unregister_edge(edge)

    cpdef void debug_inventory(self, Edge edge):
        print("Inventory on Edge:", edge.id_str)
        # ... (omitted for brevity, can be added if needed)

cdef class Region:
    def __init__(self, int weight, list dual_perimeter):
        global _region_counter
        self.weight = weight
        self.dual_perimeter = dual_perimeter
        self.dual_perimeter_edges = {edge for edge, _ in dual_perimeter}
        
        self.region_edges = set()
        self.bridge_set = set()
        self.bridge_to_region_edge_map = dict()
        
        self.id_str = f'R{_region_counter:x}'
        _region_counter += 1
        
        self.calc_point()

    cpdef void register_edge(self, Edge edge):
        self.region_edges.add(edge)
        self.bridge_set.add(edge.twin_graph_edge)
        self.bridge_to_region_edge_map[edge.twin_graph_edge] = edge

    cpdef void unregister_edge(self, Edge edge):
        self.region_edges.remove(edge)
        self.bridge_set.remove(edge.twin_graph_edge)
        del self.bridge_to_region_edge_map[edge.twin_graph_edge]

    cpdef void generate_cc_region_edges(self):
        self.cc_region_edges = []
        cdef QuadEdge edge
        for edge, _ in reversed(self.dual_perimeter):
            if edge in self.bridge_to_region_edge_map:
                self.cc_region_edges.append(self.bridge_to_region_edge_map[edge])

    cpdef void calc_point(self):
        cdef double x_sum = 0.0
        cdef double y_sum = 0.0
        cdef double scaling = 0.0
        cdef QuadEdge edge
        cdef EdgeDir dir
        cdef Vert src
        
        for edge, dir in self.dual_perimeter:
            src, _ = edge.get_dual_vert_pair(dir)
            if src.role != VertRole.DUAL_EXTERIOR:
                x_sum += src.point.x
                y_sum += src.point.y
                scaling += 1.0

        if len(self.dual_perimeter) == 0:
            self.point = Point(0.0, 0.0)
        elif scaling > 0:
            self.point = Point(x_sum / scaling, y_sum / scaling)
        else:
            self.point = Point(0.0, 0.0)

    cpdef set get_perimeter_verts(self):
        cdef set verts = set()
        cdef QuadEdge edge
        cdef EdgeDir dir
        cdef Vert src
        for edge, dir in self.dual_perimeter:
            src, _ = edge.get_dual_vert_pair(dir)
            verts.add(src)
        return verts

    cpdef tuple get_interior_lining_verts(self):
        if len(self.dual_perimeter) == 0:
            return ({e.get_dual_dest_from(self.region_tree.graph.external_dual_vert)[0] 
                    for e in self.region_tree.graph.external_dual_vert.cc_edges}, 
                    set([self.region_tree.graph.external_dual_vert]))

        cdef set full_perimeter_edges = self.dual_perimeter_edges
        cdef set lining_verts = set()
        cdef set perimeter_verts = set()
        cdef QuadEdge edge, working_edge
        cdef EdgeDir dir
        cdef Vert rotary_center

        for edge, dir in self.dual_perimeter:
            _, rotary_center = edge.get_dual_vert_pair(dir)
            perimeter_verts.add(rotary_center)

            working_edge = edge
            while True:
                working_edge, _ = working_edge.get_dual_cc_next_edge(rotary_center)
                if working_edge in full_perimeter_edges: 
                    break

                lining_verts.add(working_edge.get_dual_dest_from(rotary_center)[0])

        return lining_verts, perimeter_verts

    cpdef Vert get_uniform_interior_vert(self):
        if len(self.dual_perimeter) == 0:
            return random.choice(self.region_tree.graph.dual_verts_list)

        cdef int PERIMETER = 0
        cdef int LINING = 1
        
        cdef set lining_verts
        cdef set perimeter_verts
        lining_verts, perimeter_verts = self.get_interior_lining_verts()
        
        cdef list perimeter_table = []
        cdef Vert v
        for v in perimeter_verts:
            perimeter_table.append(BoundaryElement(v.map_tree_vert.dfs_in, PERIMETER, v))
            perimeter_table.append(BoundaryElement(v.map_tree_vert.dfs_out, PERIMETER, v))
        for v in lining_verts:
            perimeter_table.append(BoundaryElement(v.map_tree_vert.dfs_in, LINING, v))
            perimeter_table.append(BoundaryElement(v.map_tree_vert.dfs_out, LINING, v))
            
        perimeter_table.sort(key=attrgetter('dfs_index'))
        
        cdef list dual_tree_path_pairs_temp = []
        cdef BoundaryElement leading_idx, trailing_idx
        
        for leading_idx, trailing_idx in pairwise(perimeter_table):
            if leading_idx.boundary_type == LINING and trailing_idx.boundary_type == LINING:
                dual_tree_path_pairs_temp.append((leading_idx, trailing_idx))
            elif leading_idx.boundary_type == LINING and trailing_idx.boundary_type == PERIMETER:
                dual_tree_path_pairs_temp.append((leading_idx, leading_idx))
        
        if len(perimeter_table) >= 2:
            if (<BoundaryElement>perimeter_table[-2]).boundary_type == PERIMETER and (<BoundaryElement>perimeter_table[-1]).boundary_type == LINING:
                dual_tree_path_pairs_temp.append((perimeter_table[-1], perimeter_table[-1]))

        if len(dual_tree_path_pairs_temp) == 0:
            return None

        cdef list dual_tree_path_pairs = []
        cdef int double_weight = 0
        cdef tuple cursor_pair = dual_tree_path_pairs_temp[0]
        cdef tuple pair
        
        for pair in dual_tree_path_pairs_temp[1:]:
            if (<BoundaryElement>pair[0]).dfs_index == (<BoundaryElement>cursor_pair[1]).dfs_index:
                cursor_pair = (cursor_pair[0], pair[1])
            else:
                dual_tree_path_pairs.append(cursor_pair)
                leading_idx, trailing_idx = cursor_pair
                double_weight += trailing_idx.dfs_index - leading_idx.dfs_index + 1
                cursor_pair = pair
        
        dual_tree_path_pairs.append(cursor_pair)
        leading_idx, trailing_idx = cursor_pair
        double_weight += trailing_idx.dfs_index - leading_idx.dfs_index + 1
        
        if double_weight == 0:
            return None
            
        cdef int checkpoint_range_target_idx = random.randint(0, double_weight - 1)
        
        cdef int dfs_target_idx = -1
        cdef Vert entry_vert = None
        cdef int logical_range_covered = 0
        cdef int offset_into_range
        
        for leading_idx, trailing_idx in dual_tree_path_pairs:
            if logical_range_covered + (trailing_idx.dfs_index - leading_idx.dfs_index + 1) > checkpoint_range_target_idx:
                offset_into_range = checkpoint_range_target_idx - logical_range_covered
                dfs_target_idx = leading_idx.dfs_index + offset_into_range
                entry_vert = leading_idx.vert
                break
            logical_range_covered += trailing_idx.dfs_index - leading_idx.dfs_index + 1
            
        if dfs_target_idx == -1 or entry_vert is None:
            raise ValueError("Could not locate dfs target index or entry vert.")
            
        cdef MapTreeVert current_vert = entry_vert.map_tree_vert
        cdef bint contains
        cdef MapTreeVert next_vert
        
        while True:
            if current_vert.dfs_in == dfs_target_idx or current_vert.dfs_out == dfs_target_idx:
                break
            contains = (current_vert.dfs_in <= dfs_target_idx) and (dfs_target_idx <= current_vert.dfs_out)
            if contains:
                if len(current_vert.out_edges) == 0:
                    raise ValueError("Could not go down from leaf.")
                for edge in current_vert.out_edges:
                    next_vert = edge.get_dest_from(current_vert)
                    if next_vert.dfs_in <= dfs_target_idx and dfs_target_idx <= next_vert.dfs_out:
                        current_vert = next_vert
                        break
            else:
                if current_vert.in_edge is None:
                    raise ValueError("Could not go up from root.")
                current_vert = current_vert.in_edge.get_dest_from(current_vert)
                
        return current_vert.twin_vert

    cpdef Vert get_central_region_vert(self):
        cdef set visited = self.get_perimeter_verts()
        cdef set lining_verts
        lining_verts, _ = self.get_interior_lining_verts()
        cdef object queue = deque(lining_verts)
        cdef Vert dual_vert = None
        cdef QuadEdge neighbor_edge
        cdef Vert neighbor
        
        while queue:
            dual_vert = queue.popleft()
            if dual_vert in visited:
                continue
            visited.add(dual_vert)
            for neighbor_edge in dual_vert.cc_edges:
                neighbor, _ = neighbor_edge.get_dual_dest_from(dual_vert)
                if neighbor not in visited:
                    queue.append(neighbor)
        return dual_vert

    cpdef Vert get_walk_random_interior_vert(self, int steps=100):
        cdef set perimeter_verts
        cdef Vert src_vert
        cdef Vert ext
        cdef QuadEdge edge

        if len(self.dual_perimeter) == 0:
            perimeter_verts = {self.region_tree.graph.external_dual_vert}
            ext = self.region_tree.graph.external_dual_vert
            if len(ext.cc_edges) == 0:
                return None
            edge = random.choice(ext.cc_edges)
            src_vert, _ = edge.get_dual_dest_from(ext)
        else:
            perimeter_verts = self.get_perimeter_verts()
            src_vert = self._find_interior_start_vert()
        
        if src_vert is None:
            return None
            
        return self._walk_interior_optimistic(src_vert, steps, perimeter_verts)

    cdef Vert _find_interior_start_vert(self):
        cdef int start_idx = random.randint(0, len(self.dual_perimeter) - 1)
        cdef int offset = 0
        cdef int n = len(self.dual_perimeter)
        cdef QuadEdge boundary_start_edge, boundary_next_edge, next_edge
        cdef EdgeDir boundary_start_dir, boundary_next_dir
        cdef Vert boundary_next_vert, rotary_center, candidate_interior_vert
        
        while offset < n:
            boundary_start_edge, boundary_start_dir = self.dual_perimeter[(start_idx + offset) % n]
            boundary_next_edge, boundary_next_dir = self.dual_perimeter[(start_idx + offset + 1) % n]
            _, boundary_next_vert = boundary_next_edge.get_dual_vert_pair(boundary_next_dir)
            _, rotary_center = boundary_start_edge.get_dual_vert_pair(boundary_start_dir)
            
            if rotary_center is None:
                raise ValueError("Edge with undefined dual vert.")

            while True:
                next_edge, _ = boundary_start_edge.get_dual_cc_next_edge(rotary_center)
                candidate_interior_vert, _ = next_edge.get_dual_dest_from(rotary_center)

                if candidate_interior_vert is not boundary_next_vert:
                    return candidate_interior_vert

                if next_edge == boundary_next_edge:
                    break
            offset += 1
        return None

    cdef Vert _walk_interior_optimistic(self, Vert start, int step_cap, set perimeter_verts):
        cdef Vert current = start
        cdef list edges
        cdef int degree
        cdef QuadEdge edge
        cdef Vert candidate
        cdef int i
        
        for i in range(max(0, step_cap)):
            edges = current.cc_edges
            degree = len(edges)
            if degree == 0:
                break
            edge = random.choice(edges)
            candidate, _ = edge.get_dual_dest_from(current)
            if candidate not in perimeter_verts:
                current = candidate
        return current

    cpdef bint check_center(self):
        cdef bint is_center = True
        cdef Edge edge
        cdef object opposite_weight
        cdef double midpoint_weight = self.region_tree.total_weight / 2.0
        cdef double epsilon_buffer = self.region_tree.epsilon * midpoint_weight
        
        for edge in self.region_edges:
            if edge.ab_weight_differential is None:
                raise ValueError("Edge with undefined ab_weight_differential.")
            
            opposite_weight = edge.get_opposite_side_weight(self)
            if opposite_weight is None:
                raise ValueError("Edge with undefined opposite side weight.")
            elif <int>opposite_weight >= midpoint_weight + epsilon_buffer:
                is_center = False
        
        if is_center:
            if self.validate_split_possible():
                return True
        return False

    cpdef bint validate_split_possible(self):
        if len(self.cc_region_edges) == 0:
            return True

        cdef double src_graph_weight = self.region_tree.total_weight
        cdef double midpoint_weight = src_graph_weight / 2.0
        cdef double epsilon_buffer = self.region_tree.epsilon * midpoint_weight

        cdef int trailing_idx = 0
        cdef int leading_idx = 1
        cdef object a_1_weight_obj = self.cc_region_edges[0].get_opposite_side_weight(self)
        if a_1_weight_obj is None:
            raise ValueError("Edge with undefined opposite side weight.")
        cdef int a_1_weight = <int>a_1_weight_obj
        cdef int a_2_weight = <int>src_graph_weight - self.weight - a_1_weight
        
        cdef int dropped_region_weight, added_region_weight
        cdef int n = len(self.cc_region_edges)
        
        while True:
            if a_1_weight <= midpoint_weight + epsilon_buffer and a_2_weight <= midpoint_weight + epsilon_buffer:
                return True
            
            if trailing_idx == n or (leading_idx > n and leading_idx % n == trailing_idx):
                return False
            
            if a_1_weight > midpoint_weight + epsilon_buffer:
                dropped_region_weight = <int>self.cc_region_edges[trailing_idx % n].get_opposite_side_weight(self)
                a_1_weight -= dropped_region_weight
                a_2_weight += dropped_region_weight
                trailing_idx += 1
            else:
                added_region_weight = <int>self.cc_region_edges[leading_idx % n].get_opposite_side_weight(self)
                a_1_weight += added_region_weight
                a_2_weight -= added_region_weight
                leading_idx += 1

cdef class Edge:
    def __init__(self):
        global _edge_counter
        self.id_str = f'RE{_edge_counter:x}'
        _edge_counter += 1
        self.ab_weight_differential = None

    cpdef tuple get_dest_from(self, Region src):
        if src == self.end_A:
            return (self.end_B, EdgeDir.AB)
        elif src == self.end_B:
            return (self.end_A, EdgeDir.BA)
        else:
            raise KeyError("Src for get_dest_from is not on edge.")

    cpdef object get_opposite_side_weight(self, Region src):
        if self.ab_weight_differential is None:
            return None
        
        cdef int diff = <int>self.ab_weight_differential
        cdef int total = self.region_tree.graph.total_primal_weight
        
        cdef int total_weight_a = (diff + total) // 2
        cdef int total_weight_b = (total - diff) // 2

        if src == self.end_A:
            return total_weight_b
        elif src == self.end_B:
            return total_weight_a
        else:
            raise KeyError("Src for get_opposite_side_weight is not on edge.")

    cpdef void calculate_weight_differential_from_dir(self, EdgeDir dir):
        cdef int side_A_weight, side_B_weight
        cdef Edge edge
        cdef object weight
        
        if dir == EdgeDir.BA:
            side_A_weight = self.end_A.weight
            for edge in self.end_A.region_edges:
                if edge is not self:
                    weight = edge.get_opposite_side_weight(self.end_A)
                    if weight is not None:
                        side_A_weight += <int>weight
                    else:
                        raise ValueError("Undefined weight.")
            side_B_weight = self.region_tree.graph.total_primal_weight - side_A_weight
        elif dir == EdgeDir.AB:
            side_B_weight = self.end_B.weight
            for edge in self.end_B.region_edges:
                if edge is not self:
                    weight = edge.get_opposite_side_weight(self.end_B)
                    if weight is not None:
                        side_B_weight += <int>weight
                    else:
                        raise ValueError("Undefined weight.")
            side_A_weight = self.region_tree.graph.total_primal_weight - side_B_weight

        self.ab_weight_differential = side_A_weight - side_B_weight

        if abs(<int>self.ab_weight_differential) <= self.region_tree.epsilon * self.region_tree.total_weight:
            self.region_tree.edge_centers.add(self)
