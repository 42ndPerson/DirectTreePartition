from __future__ import annotations
from typing import List, Tuple, Set, Dict, Optional, Generator
from enum import Enum
import sys

from collections import deque

from Euclid import *
import Euclid

class TwinGraph:
    __slots__ = (
        'primalVerts', 
        'total_primal_weight', 
        'dualVerts', 
        'external_dual_vert',
        'edges',
        'lowerXBound', 
        'upperXBound', 
        'lowerYBound', 
        'upperYBound',
        'primalMap',
        'dualMap',
        'animating',
        'animation_tracks'
    )

    primalVerts: Set[TwinGraph.Vert]
    total_primal_weight: int
    dualVerts: Set[TwinGraph.Vert]
    external_dual_vert: TwinGraph.Vert
    edges: Set[TwinGraph.QuadEdge]

    lowerXBound: float
    upperXBound: float
    lowerYBound: float
    upperYBound: float

    primalMap: MapTree
    dualMap: MapTree

    animating: bool
    animation_tracks: List[List[List[Tuple[TwinGraph.QuadEdge, TwinGraph.VertRole, TwinGraph.EdgeDir, int]]]]

    def __init__(self, points: List[Point], weights: List[int], edgeIdxs: List[Tuple[int,int]]) -> None:
        assert len(points) == len(weights)
        zippedPointWeights = zip(points, weights)

        # Create key data structures
        ordered_primal_verts = [TwinGraph.Vert(point, weight, TwinGraph.VertRole.PRIMAL) for point, weight in zippedPointWeights]

        self.primalVerts = set(ordered_primal_verts)
        self.dualVerts = set()
        self.edges = {TwinGraph.QuadEdge(ordered_primal_verts[idxA], ordered_primal_verts[idxB]) for idxA, idxB in edgeIdxs} # Link verts at index pairs from edgeIdxs

        # Note bounds
        self.lowerXBound = min(points, key=lambda a: a.x).x
        self.upperXBound = max(points, key=lambda a: a.x).x
        self.lowerYBound = min(points, key=lambda a: a.y).y
        self.upperYBound = max(points, key=lambda a: a.y).y

        # Register edges with vertices and collect total weight
        self.total_primal_weight = 0
        vertEdgeDict = {}
        for edge in self.edges:
            vertA = edge.primal_A
            vertB = edge.primal_B
            if vertA not in vertEdgeDict:
                vertEdgeDict[vertA] = []
            if vertB not in vertEdgeDict:
                vertEdgeDict[vertB] = []
            vertEdgeDict[vertA].append(edge)
            vertEdgeDict[vertB].append(edge)
        for vert, vertEdges in vertEdgeDict.items():
            vert.register_edges(vertEdges)
            self.total_primal_weight += vert.weight

        # Aux Data
        # Animation
        self.animating = True
        self.animation_tracks = [[[]],[[]],[[]],[[]]] # [0] = dual construction, [1] = primal map tree, [2] = dual map tree, [3] = dual edge annotation

        # Construct dual graph
        self.construct_dual()

        # Construct map trees
        sys.setrecursionlimit(2*len(self.primalVerts)) # Increase recursion limit to allow deep recursions in large graphs
        self.primalMap = self.generate_map_tree(next(iter(self.primalVerts)), TwinGraph.VertRole.PRIMAL) #TODO: Find better start vert for primal
        self.dualMap = self.generate_map_tree(next(iter(self.dualVerts)), TwinGraph.VertRole.DUAL)

        # Annotate dual edges for primal counting
        self.annotate_dual_edges_for_primal_counting()

        # Populate retrieval caches
        # for edge in self.edges: # TAG: Profiling
        #     vert1, vert2 = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
        #     edge.dual_cc_next_retrieval_cache = {
        #         id(vert1): edge.get_dual_cc_next_edge(vert1),
        #         id(vert2): edge.get_dual_cc_next_edge(vert2)
        #     }

            
    def construct_dual(self) -> None:
        edges_AB: Set[TwinGraph.QuadEdge] = set(self.edges)
        edges_BA: Set[TwinGraph.QuadEdge] = set(self.edges)

        while len(edges_AB) > 0:
            src_edge = next(iter(edges_AB))
            self.construct_face_loop(
                src_edge, 
                TwinGraph.EdgeDir.AB, 
                edges_AB, 
                edges_BA
            )
        while len(edges_BA) > 0:
            src_edge = next(iter(edges_BA))
            self.construct_face_loop(
                src_edge, 
                TwinGraph.EdgeDir.BA, 
                edges_AB, 
                edges_BA
            )
        # TODO: Confirm generation of self-edges
            
        # NOTE: Adds multiplicative factore O(n d logd) where d = max degree
        for dual_vert in self.dualVerts:
            dual_vert.register_edges([]) # Triggers re-sort and link

    def construct_face_loop(
            self,
            src_edge: TwinGraph.QuadEdge,
            src_dir: TwinGraph.EdgeDir,
            edges_AB: Set[TwinGraph.QuadEdge],
            edges_BA: Set[TwinGraph.QuadEdge]):

        # Dual Vert
        dual_vert = TwinGraph.Vert(Point(0,0), 0, TwinGraph.VertRole.DUAL)

        # Info for new vert
        loop_edges: List[TwinGraph.QuadEdge] = []
        x_centroid_avg: float = 0
        y_centroid_avg: float = 0
        total_angle: float = 0
        
        # Record vertex for face loop completion and load starting vert and edge (with dir)
        root_vert, current_vert = src_edge.get_primal_vert_pair(src_dir)
        current_edge: TwinGraph.QuadEdge = src_edge
        current_dir: TwinGraph.EdgeDir = src_dir

        # Process edges
        idx = 1
        while True:
            # Record consumed edge
            if current_dir == TwinGraph.EdgeDir.AB:
                edges_AB.remove(current_edge)
            if current_dir == TwinGraph.EdgeDir.BA:
                edges_BA.remove(current_edge)
            loop_edges.append(current_edge)

            # Link edges to dual vert
            if current_dir == TwinGraph.EdgeDir.AB:
                current_edge.dual_AB = dual_vert
            if current_dir == TwinGraph.EdgeDir.BA:
                current_edge.dual_BA = dual_vert

            # Update centroid averages
            x_centroid_avg = x_centroid_avg + (current_vert.point.x - x_centroid_avg) / idx
            y_centroid_avg = y_centroid_avg + (current_vert.point.y - y_centroid_avg) / idx

            # Append list of all edges so far processed
            if self.animating:
                # Append a list of (QuadEdge, VertRole, EdgeDir) tuples for animation tracking
                edge_tuples = []
                for edge in self.edges.difference(edges_AB) | self.edges.difference(edges_BA):
                    # You may want to adjust VertRole and EdgeDir as appropriate for your animation logic
                    edge_tuples.append((edge, TwinGraph.VertRole.PRIMAL, TwinGraph.EdgeDir.AB, 0))
                self.animation_tracks[0].append(edge_tuples)

            # Break if loop is complete
            if current_vert == root_vert:
                # TODO: Check if terminates before checking last vert
                break

            # Get next edge and vert and update total angle
            prev_angle = current_edge.get_primal_rad_along(current_dir)
            current_edge, _ = current_edge.get_primal_cc_next_edge(current_vert) # Gets next edge proceeding counter-clockwise
            current_vert, current_dir = current_edge.get_primal_dest_from(current_vert)

            total_angle += Euclid.Point.normalized_angle_diff(current_edge.get_primal_rad_along(current_dir), prev_angle)

            # Idx Update
            idx += 1
        
        # Set Dual Vert Centroid and Edges
        dual_vert.point = Point(x_centroid_avg, y_centroid_avg)
        dual_vert.cc_edges = loop_edges # Dual edges do not necessarily have their other end at thsi point, so they cannot be sorted until the dual is built
        # TODO: Validate reason for getting 0
        if total_angle <= 0: # Detects counter-clockwise looping at edge of graph compared to clockwise interior looping
            dual_vert.role = TwinGraph.VertRole.DUAL_EXTERIOR
            self.external_dual_vert = dual_vert
        self.dualVerts.add(dual_vert)

    def generate_map_tree(self, root_vert: TwinGraph.Vert, role: TwinGraph.VertRole) -> TwinGraph.MapTree:
        """
        Generates a spanning tree (MapTree) over either the primal or dual graph using DFS.
        :param root_vert: The root vertex for the tree (must be in the selected role set).
        :param role: TwinGraph.VertRole.PRIMAL or TwinGraph.VertRole.DUAL
        :return: TwinGraph.MapTree instance
        """

        if role.is_primal():
            verts = self.primalVerts
            get_edges = lambda v: v.cc_edges
            get_dest_dir = lambda e, v: e.get_primal_dest_from(v)
        elif role.is_dual():
            verts = self.dualVerts
            get_edges = lambda v: v.cc_edges
            get_dest_dir = lambda e, v: e.get_dual_dest_from(v)
        else:
            raise ValueError("Invalid role for map tree generation.")

        tree_verts: Set[TwinGraph.Vert] = set(verts)
        visited_verts: Set[TwinGraph.Vert] = set()
        visited_edges: Set[TwinGraph.QuadEdge] = set()
        root: TwinGraph.MapTree.MapTreeVert = root_vert.map_tree_vert

        def dfs(
                parent_vert: TwinGraph.Vert, 
                parent_tree_vert: TwinGraph.MapTree.MapTreeVert, 
                index: int) -> int: # Returns next available index
            visited_verts.add(parent_vert)
            parent_tree_vert.dfs_in = index
            index += 1

            for edge in get_edges(parent_vert):
                dest: TwinGraph.Vert
                edge_dir: TwinGraph.EdgeDir 
                dest, edge_dir = get_dest_dir(edge, parent_vert)

                if dest not in visited_verts and dest in tree_verts:
                    child_tree_vert: TwinGraph.MapTree.MapTreeVert = dest.map_tree_vert
                    map_tree_edge = parent_tree_vert.add_child(child_tree_vert, edge, edge_dir)
                    edge.map_tree_edge = map_tree_edge

                    visited_edges.add(edge)
                    if self.animating:
                        self.animation_tracks[1 if role == TwinGraph.VertRole.PRIMAL else 2].append([(visited_edge, role, edge_dir, 0) for visited_edge in visited_edges])

                    index = dfs(dest, child_tree_vert, index)

            parent_tree_vert.dfs_out = index
            return index + 1

        dfs(root_vert, root, 0)
        return TwinGraph.MapTree(root)
    
    def map_tree_reverse_dfs_iterator(self, tree: TwinGraph.MapTree) -> Generator[TwinGraph.MapTree.MapTreeVert, None, None]:
        def dfs(vert: TwinGraph.MapTree.MapTreeVert):
            for edge in vert.out_edges:
                yield from dfs(edge.child)
            yield vert

        yield from dfs(tree.root)

    def annotate_dual_edges_for_primal_counting(self) -> None:
        if not hasattr(self, 'external_dual_vert'):
            raise ValueError("Dual graph not constructed, cannot annotate dual edges.")

        for vert in self.map_tree_reverse_dfs_iterator(self.primalMap):
            # Track sum of out edges and direction of in edge
            cycle_sum = 0
            in_edge_dir = None

            # Initialize src_vert as the dual vert on the cc edge pair to get clockwise rotation
            #  The system could work as well with counter-clockwise, but I have used clockwise elsewhere
            #  and the convention seems useful
            src_vert = vert.twin_vert.cc_edges[-1].get_common_dual_vert(vert.twin_vert.cc_edges[0])
            if src_vert is None:
                raise ValueError("Could not find common dual vert for cc edge pair.")
            # In special case of of 2-sided face, cc_edges contains inadequaite info to determine a 
            #  counter-clockwise rotation, so we check manually
            if len(vert.twin_vert.cc_edges) == 2:
                incoming_edge = vert.twin_vert.cc_edges[0] # Select 0 for incoming to reverse cc to clockwise
                outgoing_edge = vert.twin_vert.cc_edges[-1] # Select -1 for outgoing to reverse cc to clockwise

                if outgoing_edge != incoming_edge.get_dual_cc_next_edge(src_vert)[0]:
                    # If outgoing edge is not the cc next edge from src_vert along incoming_edge, we have counter-clockwise rotation
                    src_vert, _ = incoming_edge.get_dual_dest_from(src_vert) # Swap to get clockwise rotation

            for edge in reversed(vert.twin_vert.cc_edges): # Reversed to get clockwise rotation
                src_vert, dir = edge.get_dual_dest_from(src_vert)

                if self.animating:
                    # Append a list of (QuadEdge, VertRole, EdgeDir) tuples for animation tracking
                    self.animation_tracks[3].append([(edge, TwinGraph.VertRole.DUAL, TwinGraph.EdgeDir.AB, 1 if dir == TwinGraph.EdgeDir.AB else 0)])

                if edge.dual_AB_annotation is None:
                    in_edge_twin = vert.in_edge.twin_edge if vert.in_edge is not None else None
                    if edge != in_edge_twin: # All out edges will have already been annotated due to reverse_dfs order
                        edge.dual_AB_annotation = 0
                    else:
                        in_edge_dir = dir
                else:
                    # Maintain sign of cycle edges
                    if dir == TwinGraph.EdgeDir.AB:
                        cycle_sum += edge.dual_AB_annotation
                    if dir == TwinGraph.EdgeDir.BA:
                        cycle_sum -= edge.dual_AB_annotation

            # Set the in edge to balance the cycle sum to 1 (the number of primal verts in the face)
            if vert.in_edge is not None:
                if in_edge_dir is None:
                    raise ValueError("In edge direction not found for dual edge annotation.")
                if in_edge_dir == TwinGraph.EdgeDir.AB:
                    vert.in_edge.twin_edge.dual_AB_annotation = vert.twin_vert.weight - cycle_sum # TODO: Confirm the weighted version still works as intended
                if in_edge_dir == TwinGraph.EdgeDir.BA:
                    vert.in_edge.twin_edge.dual_AB_annotation = -vert.twin_vert.weight + cycle_sum

    def count_primal_verts_within_perim(self, perim: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]) -> int:
        """
        Counts the total weight of primal vertices enclosed by a given perimeter of dual edges.
        Assumes the perimeter is oriented clockwise.
        Treats an empty perimeter as enclosing the whole graph.
        :param perim: List of (QuadEdge, EdgeDir) tuples defining the perimeter. (List must 
        not be de-dupped to allow pendants in perim to self-cancel)
        """

        # Empty perimeter denotes the whole graph
        if len(perim) == 0:
            return self.total_primal_weight
        
        perim_weight = 0
        for edge, dir in perim:
            if edge.dual_AB_annotation is None:
                raise ValueError("Dual edge not annotated for primal counting.")
            if dir == TwinGraph.EdgeDir.AB:
                perim_weight += edge.dual_AB_annotation
            if dir == TwinGraph.EdgeDir.BA:
                perim_weight -= edge.dual_AB_annotation

        if perim_weight < 0:
            return perim_weight%self.total_primal_weight
        elif perim_weight == 0:
            # This case accounts for perimeters that enclose no area, such as pendants that self-cancel
            #  This case is somewhat fraught however, as it is ambiguous whether the perimeter encloses all or none of the graph
            #  This should not be a problem for this project, as we should never clockwise enclose no area, but it could cause confusion in other applications
            return self.total_primal_weight
        elif perim_weight <= self.total_primal_weight:
            return perim_weight
        else:
            raise ValueError("Perimeter weight exceeds total primal weight, something went wrong.")

    def get_verts_within_radius(self, src: Point, radius: float, role: TwinGraph.VertRole) -> List[TwinGraph.Vert]:
        matches = []

        # TODO: Integrate spatial lookup
        if role.is_primal():
            for vert in self.primalVerts:
                if Point.dist(src, vert.point) < radius:
                    matches.append(vert)
        if role.is_dual():
            for vert in self.dualVerts:
                if Point.dist(src, vert.point) < radius:
                    matches.append(vert)

        return matches

    def get_closest_vert(self, src: Point, role: TwinGraph.VertRole) -> Optional[TwinGraph.Vert]:
        max_radius = int(max(self.upperXBound - self.lowerXBound, self.upperYBound - self.lowerYBound))
        for radius in range(max_radius):
            close = self.get_verts_within_radius(src, radius, role)
            if len(close) > 0:
                closest = min(close, key=lambda vert: Point.dist(src, vert.point))
                return closest

        return None

    class EdgeDir(Enum):
        AB = 1
        BA = 2

        def reverse(self) -> TwinGraph.EdgeDir:
            if self == TwinGraph.EdgeDir.AB:
                return TwinGraph.EdgeDir.BA
            if self == TwinGraph.EdgeDir.BA:
                return TwinGraph.EdgeDir.AB

    class VertRole(Enum):
        PRIMAL = 1
        DUAL = -1
        DUAL_EXTERIOR = -2

        def is_primal(self) -> bool:
            return self.value > 0

        def is_dual(self) -> bool:
            return self.value < 0

    class Vert:
        __slots__ = (
            'point',
            'weight',
            'role',
            'cc_edges',
            'map_tree_vert',
            'id_str'
        )

        point: Point
        weight: int
        role: TwinGraph.VertRole
        cc_edges: List[TwinGraph.QuadEdge]

        map_tree_vert: TwinGraph.MapTree.MapTreeVert

        # Instance labeling
        index: int = 0
        id_str: str

        # Counterclockwise edges must be sorted
        def __init__(self, point: Point, weight: int, role: TwinGraph.VertRole, counterclockwiseEdges: List[TwinGraph.QuadEdge]=[]) -> None:
            # TODO: Add assert for edge sorting
            self.point = point
            self.weight = weight
            self.role = role
            self.cc_edges = counterclockwiseEdges
            self.link_edges()

            self.map_tree_vert = TwinGraph.MapTree.MapTreeVert(self)

            self.id_str = ("P" if role.is_primal() else "D") + f'{TwinGraph.Vert.index:x}' # Hexadecimal ID for easier reading in debug
            TwinGraph.Vert.index += 1

        # Edges must already have links to calling vertex
        def register_edges(self, edges: List[TwinGraph.QuadEdge]) -> None:
            # TODO: Add assert for edge linking
            edges = self.cc_edges + edges # TODO: Sorted insertion would be faster than re-sorting all
            if self.role.is_primal():
                edges.sort(key=lambda edge: edge.get_primal_rad_from(self)[0])
            if self.role.is_dual():
                edges.sort(key=lambda edge: edge.get_dual_rad_from(self)[0]) # TODO: Switch to system pulling order from face loops, as this solution relies on dual vert point accuracy

                # Exterior dual vert must have edges in reverse order
                if self.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                    edges.reverse() 

            self.cc_edges = edges
            self.link_edges()

        # Edges must already have links to calling vertex
        # This should already be a given, but self.cc_edges must be sorted
        def link_edges(self):
            # TODO: Add assert for edge linking
            list_len = len(self.cc_edges)
            for i, edge in enumerate(self.cc_edges):
                if self.role.is_primal():
                    if self == edge.primal_A:
                        self.cc_edges[i%list_len].primal_A_cc_next = self.cc_edges[(i+1)%list_len]
                        self.cc_edges[i%list_len].primal_A_cc_prev = self.cc_edges[(i-1)%list_len]
                    if self == edge.primal_B:
                        self.cc_edges[i%list_len].primal_B_cc_next = self.cc_edges[(i+1)%list_len]
                        self.cc_edges[i%list_len].primal_B_cc_prev = self.cc_edges[(i-1)%list_len]
                if self.role.is_dual():
                    if self == edge.dual_AB:
                        self.cc_edges[i%list_len].dual_AB_cc_next = self.cc_edges[(i+1)%list_len]
                        self.cc_edges[i%list_len].dual_AB_cc_prev = self.cc_edges[(i-1)%list_len]
                    if self == edge.dual_BA:
                        self.cc_edges[i%list_len].dual_BA_cc_next = self.cc_edges[(i+1)%list_len]
                        self.cc_edges[i%list_len].dual_BA_cc_prev = self.cc_edges[(i-1)%list_len]

    # QuadEdge contains both primal and dual edge info
    class QuadEdge:
        __slots__ = (
            'primal_A',
            'primal_B',
            'primal_A_cc_next',
            'primal_A_cc_prev',
            'primal_B_cc_next',
            'primal_B_cc_prev',
            'dual_AB',
            'dual_BA',
            'dual_AB_cc_next',
            'dual_AB_cc_prev',
            'dual_BA_cc_next',
            'dual_BA_cc_prev',
            'dual_AB_annotation',
            'map_tree_edge',
            'id_str',
            # 'dual_cc_next_retrieval_cache', # TAG: Profiling
        )

        primal_A: TwinGraph.Vert
        primal_B: TwinGraph.Vert
        primal_A_cc_next: TwinGraph.QuadEdge
        primal_A_cc_prev: TwinGraph.QuadEdge
        primal_B_cc_next: TwinGraph.QuadEdge
        primal_B_cc_prev: TwinGraph.QuadEdge

        dual_AB: TwinGraph.Vert
        dual_BA: TwinGraph.Vert
        dual_AB_cc_next: TwinGraph.QuadEdge
        dual_AB_cc_prev: TwinGraph.QuadEdge
        dual_BA_cc_next: TwinGraph.QuadEdge
        dual_BA_cc_prev: TwinGraph.QuadEdge

        dual_AB_annotation: Optional[int] # Annotation for primal vert counting
        map_tree_edge: TwinGraph.MapTree.MapTreeEdge

        # dual_cc_next_retrieval_cache: Dict[int, Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] # TAG: Profiling

        # Instance labeling
        index: int = 0
        id_str: str

        def __init__(
                self,
                primal_A: TwinGraph.Vert,
                primal_B: TwinGraph.Vert
            ) -> None:
            self.primal_A = primal_A
            self.primal_B = primal_B
            # dual_AB and dual_BA set during dual network construction at graph level
            # self.dual_AB = None
            # self.dual_BA = None

            self.dual_AB_annotation = None # Indicates unannotated

            self.id_str = f'E{TwinGraph.QuadEdge.index:x}' # Hexadecimal ID for easier reading in debug
            TwinGraph.QuadEdge.index += 1

        # Get the ordered ends of an edge based on an edge direction
        def get_primal_vert_pair(self, dir: TwinGraph.EdgeDir) -> Tuple[TwinGraph.Vert, TwinGraph.Vert]:
            if dir == TwinGraph.EdgeDir.AB:
                return (self.primal_A, self.primal_B)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.primal_B, self.primal_A)
        def get_dual_vert_pair(self, dir: TwinGraph.EdgeDir) -> Tuple[TwinGraph.Vert, TwinGraph.Vert]:
            if dir == TwinGraph.EdgeDir.AB:
                return (self.dual_AB, self.dual_BA)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.dual_BA, self.dual_AB)

        # Get the vertex on the other end of an edge along with directional annotation
        def get_primal_dest_from(self, src: TwinGraph.Vert) -> Tuple[TwinGraph.Vert, TwinGraph.EdgeDir]:
            if src is self.primal_A:
                return (self.primal_B, TwinGraph.EdgeDir.AB)
            elif src is self.primal_B:
                return (self.primal_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_primal_dest_from is not on edge.")
        def get_dual_dest_from(self, src: TwinGraph.Vert) -> Tuple[TwinGraph.Vert, TwinGraph.EdgeDir]:
            if src is self.dual_AB:
                return (self.dual_BA, TwinGraph.EdgeDir.AB)
            elif src is self.dual_BA:
                return (self.dual_AB, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_dual_dest_from is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def get_primal_rad_from(self, src: TwinGraph.Vert) -> Tuple[float, TwinGraph.EdgeDir]:
            dest, dir = self.get_primal_dest_from(src)
            angle = Point.src_dest_rad(src.point, dest.point)

            angle %= 2 * math.pi # Normalize angle to [0, 2pi]
            return (angle, dir)
        def get_dual_rad_from(self, src: TwinGraph.Vert) -> Tuple[float, TwinGraph.EdgeDir]: # TODO: Confirm no use in general project, as it relies on dual vert point accuracy
            dest, dir = self.get_dual_dest_from(src)
            if src is dest:
                raise ValueError("Src and dest are the same in get_dual_rad_from, cannot compute angle.")
            
            angle = Point.src_dest_rad(src.point, dest.point)
            if src.role == TwinGraph.VertRole.DUAL_EXTERIOR or dest.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                angle += math.pi # Flip angle if one of the verts is the exterior vert

                # Offset angle to land between primal verts to avoid edge overlap
                # TODO: Develop more robust solution
                primalA, primalB = self.get_primal_vert_pair(dir)
                primal_average = Point(
                    (primalA.point.x + primalB.point.x) / 2,
                    (primalA.point.y + primalB.point.y) / 2
                )
                angle_to_primal_average = Point.src_dest_rad(src.point if src.role == TwinGraph.VertRole.DUAL_EXTERIOR else dest.point, primal_average)
                angle_diff = angle_to_primal_average - angle
                angle += angle_diff # Offset by 1/5 of the angle difference

            angle %= 2 * math.pi # Normalize angle to [0, 2pi]
            return (angle, dir)

        # Get the angle of the edge traveling in a given dir  
        def get_primal_rad_along(self, dir: TwinGraph.EdgeDir) -> float:
            return self.get_primal_rad_from(self.primal_A if dir == TwinGraph.EdgeDir.AB else self.primal_B)[0]
        def get_dual_rad_along(self, dir: TwinGraph.EdgeDir) -> float: # TODO: Confirm no use in general project, as it relies on dual vert point accuracy
            return self.get_dual_rad_from(self.dual_AB if dir == TwinGraph.EdgeDir.AB else self.dual_BA)[0]
        
        # Get the cc next edge from vert
        def get_primal_cc_next_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert is self.primal_A:
                _, dir = self.primal_A_cc_next.get_primal_dest_from(vert)
                return (self.primal_A_cc_next, dir)
            if vert is self.primal_B:
                _, dir = self.primal_B_cc_next.get_primal_dest_from(vert)
                return (self.primal_B_cc_next, dir)
            raise KeyError("Vert for get_primal_cc_next_edge is not on edge.")
        def get_dual_cc_next_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert is self.dual_AB:
                _, dir = self.dual_AB_cc_next.get_dual_dest_from(vert)
                return (self.dual_AB_cc_next, dir)
            if vert is self.dual_BA:
                _, dir = self.dual_BA_cc_next.get_dual_dest_from(vert)
                return (self.dual_BA_cc_next, dir)
            raise KeyError("Vert for get_dual_cc_next_edge is not on edge.")
            
        # Get the cc prev edge from vert
        def get_primal_cc_prev_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert is self.primal_A:
                _, dir = self.primal_A_cc_prev.get_primal_dest_from(vert)
                return (self.primal_A_cc_prev, dir)
            if vert is self.primal_B:
                _, dir = self.primal_B_cc_prev.get_primal_dest_from(vert)
                return (self.primal_B_cc_prev, dir)
            raise KeyError("Vert for get_primal_cc_prev_edge is not on edge.")
        def get_dual_cc_prev_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert is self.dual_AB:
                _, dir = self.dual_AB_cc_prev.get_dual_dest_from(vert)
                return (self.dual_AB_cc_prev, dir)
            if vert is self.dual_BA:
                _, dir = self.dual_BA_cc_prev.get_dual_dest_from(vert)
                return (self.dual_BA_cc_prev, dir)
            raise KeyError("Vert for get_dual_cc_prev_edge is not on edge.")
        
        # Get common vert of two edges, if it exists
        def get_common_primal_vert(self, other: TwinGraph.QuadEdge) -> Optional[TwinGraph.Vert]:
            if self.primal_A is other.primal_A or self.primal_A is other.primal_B:
                return self.primal_A
            if self.primal_B is other.primal_A or self.primal_B is other.primal_B:
                return self.primal_B
            return None
        def get_common_dual_vert(self, other: TwinGraph.QuadEdge) -> Optional[TwinGraph.Vert]:
            if self.dual_AB is other.dual_AB or self.dual_AB is other.dual_BA:
                return self.dual_AB
            if self.dual_BA is other.dual_AB or self.dual_BA is other.dual_BA:
                return self.dual_BA
            return None

    class MapTree:
        __slots__ = (
            'root',
            'verts',
            'edges',
            'quad_edges'
        )

        root: TwinGraph.MapTree.MapTreeVert
        verts: Set[TwinGraph.MapTree.MapTreeVert]
        edges: Set[TwinGraph.MapTree.MapTreeEdge]
        quad_edges: Set[TwinGraph.QuadEdge]

        def __init__(self, root: MapTreeVert):
            self.root = root
            self.verts = set()
            self.edges = set()
            self.quad_edges = set()
            self.collect_tree(root)

        def collect_tree(self, vert):
            self.verts.add(vert)
            for edge in vert.out_edges:
                self.edges.add(edge)
                if edge.child not in self.verts:
                    self.collect_tree(edge.child)

        class MapTreeVert:
            __slots__ = (
                'twin_vert',
                'dfs_in',
                'dfs_out',
                'out_edges',
                'in_edge',
            )

            twin_vert: TwinGraph.Vert
            dfs_in: int
            dfs_out: int
            out_edges: List[TwinGraph.MapTree.MapTreeEdge]
            in_edge: Optional[TwinGraph.MapTree.MapTreeEdge]

            def __init__(self, twin_vert: TwinGraph.Vert):
                self.twin_vert = twin_vert
                self.dfs_in = -1
                self.dfs_out = -1
                self.out_edges = []
                self.in_edge = None

            def add_child(self, child: TwinGraph.MapTree.MapTreeVert, twin_edge: TwinGraph.QuadEdge, edge_dir: TwinGraph.EdgeDir) -> TwinGraph.MapTree.MapTreeEdge:
                edge = TwinGraph.MapTree.MapTreeEdge(self, child, twin_edge, edge_dir)
                self.out_edges.append(edge)
                child.in_edge = edge

                return edge

        class MapTreeEdge:
            __slots__ = (
                'parent',
                'child',
                'twin_edge',
                'edge_dir',
            )

            parent: TwinGraph.MapTree.MapTreeVert
            child: TwinGraph.MapTree.MapTreeVert
            twin_edge: TwinGraph.QuadEdge
            edge_dir: TwinGraph.EdgeDir

            def __init__(
                    self, 
                    parent: TwinGraph.MapTree.MapTreeVert, 
                    child: TwinGraph.MapTree.MapTreeVert, 
                    twin_edge: TwinGraph.QuadEdge, 
                    edge_dir: TwinGraph.EdgeDir
                    ):
                self.parent = parent
                self.child = child
                self.twin_edge = twin_edge
                self.edge_dir = edge_dir

            def get_dest_from(self, src: TwinGraph.MapTree.MapTreeVert) -> TwinGraph.MapTree.MapTreeVert:
                if src is self.parent:
                    return self.child
                elif src is self.child:
                    return self.parent
                else:
                    raise KeyError("Src for get_dest_from is not on edge.")