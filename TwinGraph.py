from __future__ import annotations
from typing import List, Tuple, Set, Dict, Optional
from enum import Enum

from collections import deque

from Euclid import *

class TwinGraph:
    primalVerts: Set[TwinGraph.Vert]
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

    def __init__(self, points: List[Point], weights: List[float], edgeIdxs: List[Tuple[int,int]]) -> None:
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

        # Register edges with vertices
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

        # Aux Data
        # Animation
        self.animating = True
        self.animation_tracks = [[[]],[[]],[[]]] # [0] = dual construction, [1] = primal map tree, [2] = dual map tree

        # Construct dual graph
        self.construct_dual()

        # Construct map trees
        self.primalMap = self.generate_map_tree(next(iter(self.primalVerts)), TwinGraph.VertRole.PRIMAL) #TODO: Find better start vert for primal
        self.dual = self.generate_map_tree(next(iter(self.dualVerts)), TwinGraph.VertRole.DUAL)
            
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
                break

            # Get next edge and vert and update total angle
            prev_angle = current_edge.get_primal_rad_along(current_dir)
            current_edge, _ = current_edge.get_primal_cc_next_edge(current_vert) # Gets next edge proceeding counter-clockwise
            current_vert, current_dir = current_edge.get_primal_dest_from(current_vert)

            total_angle += current_edge.get_primal_rad_along(current_dir) - prev_angle

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
                    parent_tree_vert.add_child(child_tree_vert, edge, edge_dir)
                    visited_edges.add(edge)
                    if self.animating:
                        self.animation_tracks[1 if role == TwinGraph.VertRole.PRIMAL else 2].append([(visited_edge, role, edge_dir, 0) for visited_edge in visited_edges])

                    index = dfs(dest, child_tree_vert, index)

            parent_tree_vert.dfs_out = index
            return index + 1

        dfs(root_vert, root, 0)
        return TwinGraph.MapTree(root)

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
        point: Point
        weight: float
        role: TwinGraph.VertRole
        cc_edges: List[TwinGraph.QuadEdge]

        map_tree_vert: TwinGraph.MapTree.MapTreeVert

        # Instance labeling
        index: int = 0
        id_str: str

        # Counterclockwise edges must be sorted
        def __init__(self, point: Point, weight: float, role: TwinGraph.VertRole, counterclockwiseEdges: List[TwinGraph.QuadEdge]=[]) -> None:
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
                edges.sort(key=lambda edge: edge.get_dual_rad_from(self)[0])

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
                if self == edge.primal_A:
                    self.cc_edges[i%list_len].primal_A_cc_next = self.cc_edges[(i+1)%list_len]
                    self.cc_edges[i%list_len].primal_A_cc_prev = self.cc_edges[(i-1)%list_len]
                if self == edge.primal_B:
                    self.cc_edges[i%list_len].primal_B_cc_next = self.cc_edges[(i+1)%list_len]
                    self.cc_edges[i%list_len].primal_B_cc_prev = self.cc_edges[(i-1)%list_len]
                if self == edge.dual_AB:
                    self.cc_edges[i%list_len].dual_AB_cc_next = self.cc_edges[(i+1)%list_len]
                    self.cc_edges[i%list_len].dual_AB_cc_prev = self.cc_edges[(i-1)%list_len]
                if self == edge.dual_BA:
                    self.cc_edges[i%list_len].dual_BA_cc_next = self.cc_edges[(i+1)%list_len]
                    self.cc_edges[i%list_len].dual_BA_cc_prev = self.cc_edges[(i-1)%list_len]

    # QuadEdge contains both primal and dual edge info
    class QuadEdge:
        primal_A: TwinGraph.Vert
        primal_B: TwinGraph.Vert
        primal_AB_annotation: Optional[int]
        primal_BA_annotation: Optional[int]
        primal_A_cc_next: TwinGraph.QuadEdge
        primal_A_cc_prev: TwinGraph.QuadEdge
        primal_B_cc_next: TwinGraph.QuadEdge
        primal_B_cc_prev: TwinGraph.QuadEdge

        dual_AB: Optional[TwinGraph.Vert]
        dual_BA: Optional[TwinGraph.Vert]
        dual_AB_annotation: Optional[int]
        dual_BA_annotation: Optional[int]
        dual_AB_cc_next: TwinGraph.QuadEdge
        dual_AB_cc_prev: TwinGraph.QuadEdge
        dual_BA_cc_next: TwinGraph.QuadEdge
        dual_BA_cc_prev: TwinGraph.QuadEdge

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
            self.primal_AB_annotation = None
            self.primal_BA_annotation = None
            self.dual_AB = None
            self.dual_BA = None
            self.dual_AB_annotation = None
            self.dual_BA_annotation = None

            self.id_str = f'E{TwinGraph.QuadEdge.index:x}' # Hexadecimal ID for easier reading in debug
            TwinGraph.QuadEdge.index += 1

        # Get the ordered ends of an edge based on an edge direction
        def get_primal_vert_pair(self, dir: TwinGraph.EdgeDir) -> Tuple[TwinGraph.Vert, TwinGraph.Vert]:
            if dir == TwinGraph.EdgeDir.AB:
                return (self.primal_A, self.primal_B)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.primal_B, self.primal_A)
        def get_dual_vert_pair(self, dir: TwinGraph.EdgeDir) -> Tuple[TwinGraph.Vert, TwinGraph.Vert]:
            if self.dual_AB is None or self.dual_BA is None:
                raise ValueError("Dual edge not yet constructed.")
            if dir == TwinGraph.EdgeDir.AB:
                return (self.dual_AB, self.dual_BA)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.dual_BA, self.dual_AB)

        # Get the vertex on the other end of an edge along with directional annotation
        def get_primal_dest_from(self, src: TwinGraph.Vert) -> Tuple[TwinGraph.Vert, TwinGraph.EdgeDir]:
            if src == self.primal_A:
                return (self.primal_B, TwinGraph.EdgeDir.AB)
            elif src == self.primal_B:
                return (self.primal_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_primal_dest_from is not on edge.")
        def get_dual_dest_from(self, src: TwinGraph.Vert) -> Tuple[TwinGraph.Vert, TwinGraph.EdgeDir]:
            if src == self.dual_AB:
                if self.dual_BA is None:
                    raise ValueError("Dual edge not yet constructed.")
                return (self.dual_BA, TwinGraph.EdgeDir.AB)
            elif src == self.dual_BA:
                if self.dual_AB is None:
                    raise ValueError("Dual edge not yet constructed.")
                return (self.dual_AB, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_dual_dest_from is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def get_primal_rad_from(self, src: TwinGraph.Vert) -> Tuple[float, TwinGraph.EdgeDir]:
            dest, dir = self.get_primal_dest_from(src)
            angle = Point.src_dest_rad(src.point, dest.point)

            angle %= 2 * math.pi # Normalize angle to [0, 2pi]
            return (angle, dir)
        def get_dual_rad_from(self, src: TwinGraph.Vert) -> Tuple[float, TwinGraph.EdgeDir]:
            dest, dir = self.get_dual_dest_from(src)
            if src == dest:
                raise ValueError("Src and dest are the same in get_dual_rad_from, cannot compute angle.")
            
            angle = Point.src_dest_rad(src.point, dest.point)
            if src.role == TwinGraph.VertRole.DUAL_EXTERIOR or dest.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                angle += math.pi # Flip angle if one of the verts is the exterior vert

            angle %= 2 * math.pi # Normalize angle to [0, 2pi]
            return (angle, dir)

        # Get the angle of the edge traveling in a given dir  
        def get_primal_rad_along(self, dir: TwinGraph.EdgeDir) -> float:
            return self.get_primal_rad_from(self.primal_A if dir == TwinGraph.EdgeDir.AB else self.primal_B)[0]
            # if dir == TwinGraph.EdgeDir.AB:
            #     return Point.src_dest_rad(self.primal_A.point, self.primal_B.point)
            # if dir == TwinGraph.EdgeDir.BA:
            #     return Point.src_dest_rad(self.primal_B.point, self.primal_A.point)
        def get_dual_rad_along(self, dir: TwinGraph.EdgeDir) -> float:
            if self.dual_AB is None or self.dual_BA is None:
                raise ValueError("Dual vertices not set for this edge.")
            return self.get_dual_rad_from(self.dual_AB if dir == TwinGraph.EdgeDir.AB else self.dual_BA)[0]
            # if dir == TwinGraph.EdgeDir.AB:
            #     return Point.src_dest_rad(self.dual_AB.point, self.dual_BA.point)
            # if dir == TwinGraph.EdgeDir.BA:
            #     return Point.src_dest_rad(self.dual_BA.point, self.dual_AB.point)
        
        # Get the cc next edge from vert
        def get_primal_cc_next_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert == self.primal_A:
                _, dir = self.primal_A_cc_next.get_primal_dest_from(vert)
                return (self.primal_A_cc_next, dir)
            if vert == self.primal_B:
                _, dir = self.primal_B_cc_next.get_primal_dest_from(vert)
                return (self.primal_B_cc_next, dir)
            raise KeyError("Vert for get_primal_cc_next_edge is not on edge.")
        def get_dual_cc_next_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert == self.dual_AB:
                _, dir = self.dual_AB_cc_next.get_dual_dest_from(vert)
                return (self.dual_AB_cc_next, dir)
            if vert == self.dual_BA:
                _, dir = self.dual_BA_cc_next.get_dual_dest_from(vert)
                return (self.dual_BA_cc_next, dir)
            raise KeyError("Vert for get_dual_cc_next_edge is not on edge.")
            
        # Get the cc prev edge from vert
        def get_primal_cc_prev_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert == self.primal_A:
                _, dir = self.primal_A_cc_prev.get_primal_dest_from(vert)
                return (self.primal_A_cc_prev, dir)
            if vert == self.primal_B:
                _, dir = self.primal_B_cc_prev.get_primal_dest_from(vert)
                return (self.primal_B_cc_prev, dir)
            raise KeyError("Vert for get_primal_cc_prev_edge is not on edge.")
        def get_dual_cc_prev_edge(self, vert: TwinGraph.Vert) -> Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]:
            if vert == self.dual_AB:
                _, dir = self.dual_AB_cc_prev.get_dual_dest_from(vert)
                return (self.dual_AB_cc_prev, dir)
            if vert == self.dual_BA:
                _, dir = self.dual_BA_cc_prev.get_dual_dest_from(vert)
                return (self.dual_BA_cc_prev, dir)
            raise KeyError("Vert for get_dual_cc_prev_edge is not on edge.")
        
    class MapTree:
        root: TwinGraph.MapTree.MapTreeVert
        verts: Set[TwinGraph.MapTree.MapTreeVert]
        edges: Set[TwinGraph.MapTree.MapTreeEdge]

        def __init__(self, root: MapTreeVert):
            self.root = root
            self.verts = set()
            self.edges = set()
            self.collect_tree(root)

        def collect_tree(self, vert):
            self.verts.add(vert)
            for edge in vert.out_edges:
                self.edges.add(edge)
                if edge.child not in self.verts:
                    self.collect_tree(edge.child)

        class MapTreeVert:
            def __init__(self, twin_vert: TwinGraph.Vert):
                self.twin_vert = twin_vert
                self.dfs_in: int = -1
                self.dfs_out: int = -1
                self.out_edges: List[TwinGraph.MapTree.MapTreeEdge] = []
                self.in_edge: Optional[TwinGraph.MapTree.MapTreeEdge] = None

            def add_child(self, child: TwinGraph.MapTree.MapTreeVert, twin_edge: TwinGraph.QuadEdge, edge_dir: TwinGraph.EdgeDir):
                edge = TwinGraph.MapTree.MapTreeEdge(self, child, twin_edge, edge_dir)
                self.out_edges.append(edge)
                child.in_edge = edge

        class MapTreeEdge:
            def __init__(self, parent: TwinGraph.MapTree.MapTreeVert, child: TwinGraph.MapTree.MapTreeVert, twin_edge: TwinGraph.QuadEdge, edge_dir: TwinGraph.EdgeDir):
                self.parent = parent
                self.child = child
                self.twin_edge = twin_edge
                self.edge_dir = edge_dir