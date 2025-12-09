# cython: language_level=3
import sys
from libc.math cimport pi, sqrt, atan2

from CythonCore.Euclid cimport Point

# Global counters
cdef int _vert_counter = 0
cdef int _edge_counter = 0

cpdef EdgeDir reverse_dir(EdgeDir dir):
    if dir == EdgeDir.AB:
        return EdgeDir.BA
    if dir == EdgeDir.BA:
        return EdgeDir.AB
    raise ValueError("Invalid EdgeDir value")

cdef class Vert:
    def __init__(self, Point point, int weight, VertRole role, list counterclockwiseEdges=None):
        global _vert_counter
        self.point = point
        self.weight = weight
        self.role = role
        self.cc_edges = counterclockwiseEdges if counterclockwiseEdges is not None else []
        self.link_edges()

        self.map_tree_vert = MapTreeVert(self)

        self.index = _vert_counter
        self.id_str = ("P" if role > 0 else "D") + f'{self.index:x}'
        _vert_counter += 1

    cpdef void register_edges(self, list edges):
        cdef QuadEdge edge
        # edges = self.cc_edges + edges # Inefficient
        self.cc_edges.extend(edges)
        
        cdef list pairs = []
        cdef double angle
        
        if self.role > 0: # PRIMAL
            for edge in self.cc_edges:
                angle = (<QuadEdge>edge).get_primal_rad_from(self)[0]
                pairs.append((angle, edge))
            pairs.sort()
            self.cc_edges = [p[1] for p in pairs]
        elif self.role < 0: # DUAL or DUAL_EXTERIOR
            for edge in self.cc_edges:
                angle = (<QuadEdge>edge).get_dual_rad_from(self)[0]
                pairs.append((angle, edge))
            pairs.sort()
            self.cc_edges = [p[1] for p in pairs]
            
            if self.role == VertRole.DUAL_EXTERIOR:
                self.cc_edges.reverse()

        self.link_edges()

    cpdef void link_edges(self):
        cdef int list_len = len(self.cc_edges)
        cdef int i
        cdef QuadEdge edge
        
        for i in range(list_len):
            edge = self.cc_edges[i]
            if self.role > 0: # PRIMAL
                if self is edge.primal_A:
                    edge.primal_A_cc_next = self.cc_edges[(i+1)%list_len]
                    edge.primal_A_cc_prev = self.cc_edges[(i-1)%list_len]
                elif self is edge.primal_B:
                    edge.primal_B_cc_next = self.cc_edges[(i+1)%list_len]
                    edge.primal_B_cc_prev = self.cc_edges[(i-1)%list_len]
            elif self.role < 0: # DUAL
                if self is edge.dual_AB:
                    edge.dual_AB_cc_next = self.cc_edges[(i+1)%list_len]
                    edge.dual_AB_cc_prev = self.cc_edges[(i-1)%list_len]
                elif self is edge.dual_BA:
                    edge.dual_BA_cc_next = self.cc_edges[(i+1)%list_len]
                    edge.dual_BA_cc_prev = self.cc_edges[(i-1)%list_len]

cdef class QuadEdge:
    def __init__(self, Vert primal_A, Vert primal_B):
        global _edge_counter
        self.primal_A = primal_A
        self.primal_B = primal_B
        self.dual_AB_annotation = None
        
        self.id_str = f'E{_edge_counter:x}'
        _edge_counter += 1

    cpdef tuple get_primal_vert_pair(self, EdgeDir dir):
        if dir == EdgeDir.AB:
            return (self.primal_A, self.primal_B)
        if dir == EdgeDir.BA:
            return (self.primal_B, self.primal_A)
        return (None, None)

    cpdef tuple get_dual_vert_pair(self, EdgeDir dir):
        if dir == EdgeDir.AB:
            return (self.dual_AB, self.dual_BA)
        if dir == EdgeDir.BA:
            return (self.dual_BA, self.dual_AB)
        return (None, None)

    cpdef tuple get_primal_dest_from(self, Vert src):
        if src is self.primal_A:
            return (self.primal_B, EdgeDir.AB)
        elif src is self.primal_B:
            return (self.primal_A, EdgeDir.BA)
        else:
            raise KeyError("Src for get_primal_dest_from is not on edge.")

    cpdef tuple get_dual_dest_from(self, Vert src):
        if src is self.dual_AB:
            return (self.dual_BA, EdgeDir.AB)
        elif src is self.dual_BA:
            return (self.dual_AB, EdgeDir.BA)
        else:
            raise KeyError("Src for get_dual_dest_from is not on edge.")

    cpdef tuple get_primal_rad_from(self, Vert src):
        cdef Vert dest
        cdef EdgeDir dir
        cdef double angle
        dest, dir = self.get_primal_dest_from(src)
        angle = Point.src_dest_rad(src.point, dest.point)
        angle %= 2 * pi
        return (angle, dir)

    cpdef tuple get_dual_rad_from(self, Vert src):
        cdef Vert dest
        cdef EdgeDir dir
        cdef double angle
        cdef Vert primalA
        cdef Vert primalB
        cdef Point primal_average
        cdef double angle_to_primal_average
        cdef double angle_diff

        dest, dir = self.get_dual_dest_from(src)
        if src is dest:
            raise ValueError("Src and dest are the same in get_dual_rad_from")
        
        angle = Point.src_dest_rad(src.point, dest.point)
        if src.role == VertRole.DUAL_EXTERIOR or dest.role == VertRole.DUAL_EXTERIOR:
            angle += pi
            
            primalA, primalB = self.get_primal_vert_pair(dir)
            primal_average = Point(
                (primalA.point.x + primalB.point.x) / 2,
                (primalA.point.y + primalB.point.y) / 2
            )
            angle_to_primal_average = Point.src_dest_rad(src.point if src.role == VertRole.DUAL_EXTERIOR else dest.point, primal_average)
            angle_diff = angle_to_primal_average - angle
            angle += angle_diff

        angle %= 2 * pi
        return (angle, dir)

    cpdef double get_primal_rad_along(self, EdgeDir dir):
        return self.get_primal_rad_from(self.primal_A if dir == EdgeDir.AB else self.primal_B)[0]

    cpdef double get_dual_rad_along(self, EdgeDir dir):
        return self.get_dual_rad_from(self.dual_AB if dir == EdgeDir.AB else self.dual_BA)[0]

    cpdef tuple get_primal_cc_next_edge(self, Vert vert):
        cdef EdgeDir dir
        if vert is self.primal_A:
            _, dir = self.primal_A_cc_next.get_primal_dest_from(vert)
            return (self.primal_A_cc_next, dir)
        if vert is self.primal_B:
            _, dir = self.primal_B_cc_next.get_primal_dest_from(vert)
            return (self.primal_B_cc_next, dir)
        raise KeyError("Vert for get_primal_cc_next_edge is not on edge.")

    cpdef tuple get_dual_cc_next_edge(self, Vert vert):
        cdef EdgeDir dir
        if vert is self.dual_AB:
            _, dir = self.dual_AB_cc_next.get_dual_dest_from(vert)
            return (self.dual_AB_cc_next, dir)
        if vert is self.dual_BA:
            _, dir = self.dual_BA_cc_next.get_dual_dest_from(vert)
            return (self.dual_BA_cc_next, dir)
        raise KeyError("Vert for get_dual_cc_next_edge is not on edge.")

    cpdef tuple get_primal_cc_prev_edge(self, Vert vert):
        cdef EdgeDir dir
        if vert is self.primal_A:
            _, dir = self.primal_A_cc_prev.get_primal_dest_from(vert)
            return (self.primal_A_cc_prev, dir)
        if vert is self.primal_B:
            _, dir = self.primal_B_cc_prev.get_primal_dest_from(vert)
            return (self.primal_B_cc_prev, dir)
        raise KeyError("Vert for get_primal_cc_prev_edge is not on edge.")

    cpdef tuple get_dual_cc_prev_edge(self, Vert vert):
        cdef EdgeDir dir
        if vert is self.dual_AB:
            _, dir = self.dual_AB_cc_prev.get_dual_dest_from(vert)
            return (self.dual_AB_cc_prev, dir)
        if vert is self.dual_BA:
            _, dir = self.dual_BA_cc_prev.get_dual_dest_from(vert)
            return (self.dual_BA_cc_prev, dir)
        raise KeyError("Vert for get_dual_cc_prev_edge is not on edge.")

    cpdef Vert get_common_primal_vert(self, QuadEdge other):
        if self.primal_A is other.primal_A or self.primal_A is other.primal_B:
            return self.primal_A
        if self.primal_B is other.primal_A or self.primal_B is other.primal_B:
            return self.primal_B
        return None

    cpdef Vert get_common_dual_vert(self, QuadEdge other):
        if self.dual_AB is other.dual_AB or self.dual_AB is other.dual_BA:
            return self.dual_AB
        if self.dual_BA is other.dual_AB or self.dual_BA is other.dual_BA:
            return self.dual_BA
        return None

cdef class MapTreeVert:
    def __init__(self, Vert twin_vert):
        self.twin_vert = twin_vert
        self.dfs_in = -1
        self.dfs_out = -1
        self.out_edges = []
        self.in_edge = None

    cpdef MapTreeEdge add_child(self, MapTreeVert child, QuadEdge twin_edge, EdgeDir edge_dir):
        cdef MapTreeEdge edge = MapTreeEdge(self, child, twin_edge, edge_dir)
        self.out_edges.append(edge)
        child.in_edge = edge
        return edge

cdef class MapTreeEdge:
    def __init__(self, MapTreeVert parent, MapTreeVert child, QuadEdge twin_edge, EdgeDir edge_dir):
        self.parent = parent
        self.child = child
        self.twin_edge = twin_edge
        self.edge_dir = edge_dir

    cpdef MapTreeVert get_dest_from(self, MapTreeVert src):
        if src is self.parent:
            return self.child
        elif src is self.child:
            return self.parent
        else:
            raise KeyError("Src for get_dest_from is not on edge.")

cdef class MapTree:
    def __init__(self, MapTreeVert root):
        self.root = root
        self.verts = set()
        self.edges = set()
        self.quad_edges = set()
        self.collect_tree(root)

    cpdef void collect_tree(self, MapTreeVert vert):
        self.verts.add(vert)
        cdef MapTreeEdge edge
        for edge in vert.out_edges:
            self.edges.add(edge)
            if edge.child not in self.verts:
                self.collect_tree(edge.child)

cdef class TwinGraph:
    def __init__(self, list points, list weights, list edgeIdxs):
        assert len(points) == len(weights)
        
        cdef int i
        cdef Point p
        cdef int w
        cdef Vert v
        cdef list ordered_primal_verts = []
        
        for i in range(len(points)):
            p = points[i]
            w = weights[i]
            v = Vert(p, w, VertRole.PRIMAL)
            ordered_primal_verts.append(v)

        self.primalVerts = set(ordered_primal_verts)
        self.dualVerts = set()
        
        cdef int idxA, idxB
        cdef QuadEdge edge
        self.edges = set()
        for idxA, idxB in edgeIdxs:
            edge = QuadEdge(ordered_primal_verts[idxA], ordered_primal_verts[idxB])
            self.edges.add(edge)

        self.lowerXBound = min(points, key=lambda a: (<Point>a).x).x
        self.upperXBound = max(points, key=lambda a: (<Point>a).x).x
        self.lowerYBound = min(points, key=lambda a: (<Point>a).y).y
        self.upperYBound = max(points, key=lambda a: (<Point>a).y).y

        self.total_primal_weight = 0
        cdef dict vertEdgeDict = {}
        cdef Vert vertA, vertB
        
        for edge in self.edges:
            vertA = edge.primal_A
            vertB = edge.primal_B
            if vertA not in vertEdgeDict:
                vertEdgeDict[vertA] = []
            if vertB not in vertEdgeDict:
                vertEdgeDict[vertB] = []
            vertEdgeDict[vertA].append(edge)
            vertEdgeDict[vertB].append(edge)
            
        for v, vertEdges in vertEdgeDict.items():
            v.register_edges(vertEdges)
            self.total_primal_weight += v.weight

        self.animating = False
        self.animation_tracks = [[[]],[[]],[[]],[[]]]

        self.construct_dual()
        self.dual_verts_list = list(self.dualVerts)

        sys.setrecursionlimit(2*len(self.primalVerts))
        self.primalMap = self.generate_map_tree(next(iter(self.primalVerts)), VertRole.PRIMAL)
        self.dualMap = self.generate_map_tree(next(iter(self.dualVerts)), VertRole.DUAL)
        
        self.annotate_dual_edges_for_primal_counting()

    cpdef void construct_dual(self):
        cdef set edges_AB = set(self.edges)
        cdef set edges_BA = set(self.edges)
        cdef QuadEdge src_edge
        cdef Vert dual_vert

        while len(edges_AB) > 0:
            src_edge = next(iter(edges_AB))
            self.construct_face_loop(src_edge, EdgeDir.AB, edges_AB, edges_BA)
        while len(edges_BA) > 0:
            src_edge = next(iter(edges_BA))
            self.construct_face_loop(src_edge, EdgeDir.BA, edges_AB, edges_BA)
            
        for dual_vert in self.dualVerts:
            dual_vert.register_edges([])

    cpdef void construct_face_loop(self, QuadEdge src_edge, EdgeDir src_dir, set edges_AB, set edges_BA):
        cdef Vert dual_vert = Vert(Point(0,0), 0, VertRole.DUAL)
        cdef list loop_edges = []
        cdef double x_centroid_avg = 0
        cdef double y_centroid_avg = 0
        cdef double total_angle = 0
        
        cdef Vert root_vert
        cdef Vert current_vert
        root_vert, current_vert = src_edge.get_primal_vert_pair(src_dir)
        
        cdef QuadEdge current_edge = src_edge
        cdef EdgeDir current_dir = src_dir
        cdef int idx = 1
        cdef double prev_angle
        cdef list edge_tuples

        while True:
            if current_dir == EdgeDir.AB:
                edges_AB.remove(current_edge)
            if current_dir == EdgeDir.BA:
                edges_BA.remove(current_edge)
            loop_edges.append(current_edge)

            if current_dir == EdgeDir.AB:
                current_edge.dual_AB = dual_vert
            if current_dir == EdgeDir.BA:
                current_edge.dual_BA = dual_vert

            x_centroid_avg = x_centroid_avg + (current_vert.point.x - x_centroid_avg) / idx
            y_centroid_avg = y_centroid_avg + (current_vert.point.y - y_centroid_avg) / idx

            if self.animating:
                edge_tuples = []
                for edge in self.edges.difference(edges_AB) | self.edges.difference(edges_BA):
                    edge_tuples.append((edge, VertRole.PRIMAL, EdgeDir.AB, 0))
                self.animation_tracks[0].append(edge_tuples)

            if current_vert == root_vert:
                break

            prev_angle = current_edge.get_primal_rad_along(current_dir)
            current_edge, _ = current_edge.get_primal_cc_next_edge(current_vert)
            current_vert, current_dir = current_edge.get_primal_dest_from(current_vert)

            total_angle += Point.normalized_angle_diff(current_edge.get_primal_rad_along(current_dir), prev_angle)
            idx += 1
        
        dual_vert.point = Point(x_centroid_avg, y_centroid_avg)
        dual_vert.cc_edges = loop_edges
        
        if total_angle <= 0:
            dual_vert.role = VertRole.DUAL_EXTERIOR
            self.external_dual_vert = dual_vert
        self.dualVerts.add(dual_vert)

    cpdef MapTree generate_map_tree(self, Vert root_vert, VertRole role):
        cdef set verts
        if role > 0: # PRIMAL
            verts = self.primalVerts
        elif role < 0: # DUAL
            verts = self.dualVerts
        else:
            raise ValueError("Invalid role for map tree generation.")

        cdef set tree_verts = set(verts)
        cdef set visited_verts = set()
        cdef set visited_edges = set()
        cdef MapTreeVert root = root_vert.map_tree_vert

        self._dfs(root_vert, root, 0, visited_verts, visited_edges, tree_verts, role)
        return MapTree(root)

    cdef int _dfs(self, Vert parent_vert, MapTreeVert parent_tree_vert, int index, set visited_verts, set visited_edges, set tree_verts, VertRole role):
        visited_verts.add(parent_vert)
        parent_tree_vert.dfs_in = index
        index += 1

        cdef QuadEdge edge
        cdef Vert dest
        cdef EdgeDir edge_dir
        cdef MapTreeVert child_tree_vert
        cdef MapTreeEdge map_tree_edge
        cdef list edges

        edges = parent_vert.cc_edges

        for edge in edges:
            if role > 0:
                dest, edge_dir = edge.get_primal_dest_from(parent_vert)
            else:
                dest, edge_dir = edge.get_dual_dest_from(parent_vert)

            if dest not in visited_verts and dest in tree_verts:
                child_tree_vert = dest.map_tree_vert
                map_tree_edge = parent_tree_vert.add_child(child_tree_vert, edge, edge_dir)
                edge.map_tree_edge = map_tree_edge

                visited_edges.add(edge)
                if self.animating:
                    self.animation_tracks[1 if role == VertRole.PRIMAL else 2].append([(visited_edge, role, edge_dir, 0) for visited_edge in visited_edges])

                index = self._dfs(dest, child_tree_vert, index, visited_verts, visited_edges, tree_verts, role)

        parent_tree_vert.dfs_out = index
        return index + 1

    cpdef void populate_index_to_dfs_lookup(self):
        cdef int size = len(self.primalVerts) + len(self.dualVerts)
        self.idx_to_dfs_in = [0] * size
        self.idx_to_dfs_out = [0] * size

        cdef Vert vert
        for vert in self.primalVerts:
            self.idx_to_dfs_in[vert.index] = vert.map_tree_vert.dfs_in
            self.idx_to_dfs_out[vert.index] = vert.map_tree_vert.dfs_out

        for vert in self.dualVerts:
            self.idx_to_dfs_in[vert.index] = vert.map_tree_vert.dfs_in
            self.idx_to_dfs_out[vert.index] = vert.map_tree_vert.dfs_out

    def map_tree_reverse_dfs_iterator(self, MapTree tree):
        yield from self._reverse_dfs(tree.root)

    def _reverse_dfs(self, MapTreeVert vert):
        cdef MapTreeEdge edge
        for edge in vert.out_edges:
            yield from self._reverse_dfs(edge.child)
        yield vert

    cpdef void annotate_dual_edges_for_primal_counting(self):
        if not hasattr(self, 'external_dual_vert'):
            raise ValueError("Dual graph not constructed")

        cdef MapTreeVert vert
        cdef int cycle_sum
        cdef EdgeDir in_edge_dir
        cdef Vert src_vert
        cdef QuadEdge incoming_edge, outgoing_edge
        cdef QuadEdge edge
        cdef EdgeDir dir
        cdef QuadEdge in_edge_twin

        for vert in self.map_tree_reverse_dfs_iterator(self.primalMap):
            cycle_sum = 0
            in_edge_dir = EdgeDir.AB # Default, will be overwritten if in_edge exists
            
            src_vert = vert.twin_vert.cc_edges[-1].get_common_dual_vert(vert.twin_vert.cc_edges[0])
            if src_vert is None:
                raise ValueError("Could not find common dual vert")
            
            if len(vert.twin_vert.cc_edges) == 2:
                incoming_edge = vert.twin_vert.cc_edges[0]
                outgoing_edge = vert.twin_vert.cc_edges[-1]
                if outgoing_edge != incoming_edge.get_dual_cc_next_edge(src_vert)[0]:
                    src_vert, _ = incoming_edge.get_dual_dest_from(src_vert)

            for edge in reversed(vert.twin_vert.cc_edges):
                src_vert, dir = edge.get_dual_dest_from(src_vert)

                if self.animating:
                    self.animation_tracks[3].append([(edge, VertRole.DUAL, EdgeDir.AB, 1 if dir == EdgeDir.AB else 0)])

                if edge.dual_AB_annotation is None:
                    in_edge_twin = vert.in_edge.twin_edge if vert.in_edge is not None else None
                    if edge != in_edge_twin:
                        edge.dual_AB_annotation = 0
                    else:
                        in_edge_dir = dir
                else:
                    if dir == EdgeDir.AB:
                        cycle_sum += <int>edge.dual_AB_annotation
                    if dir == EdgeDir.BA:
                        cycle_sum -= <int>edge.dual_AB_annotation

            if vert.in_edge is not None:
                if in_edge_dir == EdgeDir.AB:
                    vert.in_edge.twin_edge.dual_AB_annotation = vert.twin_vert.weight - cycle_sum
                if in_edge_dir == EdgeDir.BA:
                    vert.in_edge.twin_edge.dual_AB_annotation = -vert.twin_vert.weight + cycle_sum

    cpdef int count_primal_verts_within_perim(self, list perim):
        if len(perim) == 0:
            return self.total_primal_weight
        
        cdef int perim_weight = 0
        cdef QuadEdge edge
        cdef EdgeDir dir
        
        for edge, dir in perim:
            if edge.dual_AB_annotation is None:
                raise ValueError("Dual edge not annotated")
            if dir == EdgeDir.AB:
                perim_weight += <int>edge.dual_AB_annotation
            if dir == EdgeDir.BA:
                perim_weight -= <int>edge.dual_AB_annotation

        if perim_weight < 0:
            return perim_weight % self.total_primal_weight
        elif perim_weight == 0:
            return self.total_primal_weight
        elif perim_weight <= self.total_primal_weight:
            return perim_weight
        else:
            raise ValueError("Perimeter weight exceeds total primal weight")

    cpdef list get_verts_within_radius(self, Point src, double radius, VertRole role):
        cdef list matches = []
        cdef Vert vert
        if role > 0:
            for vert in self.primalVerts:
                if Point.dist(src, vert.point) < radius:
                    matches.append(vert)
        if role < 0:
            for vert in self.dualVerts:
                if Point.dist(src, vert.point) < radius:
                    matches.append(vert)
        return matches

    cpdef Vert get_closest_vert(self, Point src, VertRole role):
        cdef int max_radius = int(max(self.upperXBound - self.lowerXBound, self.upperYBound - self.lowerYBound))
        cdef int radius
        cdef list close
        cdef Vert best_vert = None
        cdef double min_dist = float('inf')
        cdef double d
        cdef Vert v

        for radius in range(max_radius):
            close = self.get_verts_within_radius(src, float(radius), role)
            if len(close) > 0:
                best_vert = None
                min_dist = float('inf')
                for v in close:
                    d = Point.dist(src, v.point)
                    if d < min_dist:
                        min_dist = d
                        best_vert = v
                return best_vert
        return None
