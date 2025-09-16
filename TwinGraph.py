from __future__ import annotations

from collections import deque
import bisect

from enum import Enum
from Euclid import *

class TwinGraph:
    primalVerts: set(TwinGraph.Vert)
    dualVerts: set(TwinGraph.Vert)
    edges: set(TwinGraph.QuadEdge)

    lowerXBound: float
    upperXBound: float
    lowerYBound: float
    upperYBound: float

    spatial_lookup_res: float
    spatial_lookup_ref: [[TwinGraph.Vert]] # Indexing on x and then y

    animating: bool
    animation_track: [[(TwinGraph.QuadEdge, TwinGraph.VertRole, TwinGraph.EdgeDir)]]

    def __init__(self, points: [Point], weights: [float], edgeIdxs: [(int,int)]) -> None:
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
        self.animation_track = [[]]

        # Construct dual graph
        self.construct_dual()
            
    def construct_dual(self) -> None:
        edges_AB: set(TwinGraph.QuadEdge) = set(self.edges)
        edges_BA: set(TwinGraph.QuadEdge) = set(self.edges)

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
            edges_AB: set(TwinGraph.QuadEdge),
            edges_BA: set(TwinGraph.QuadEdge)):

        # Dual Vert
        dual_vert = TwinGraph.Vert(Point(0,0), 0, TwinGraph.VertRole.DUAL)

        # Info for new vert
        loop_edges: [TwinGraph.QuadEdge] = []
        x_centroid_avg: float = 0
        y_centroid_avg: float = 0
        
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
                self.animation_track.append(list(self.edges.difference(edges_AB) | self.edges.difference(edges_BA)))

            # Break if loop is complete
            if current_vert == root_vert:
                break

            # Get next edge and vert
            current_edge, _ = current_edge.get_primal_cc_next_edge(current_vert) # Gets next edge proceeding counter-clockwise
            current_vert, current_dir = current_edge.get_primal_dest_from(current_vert)

            # Idx Update
            idx += 1
        
        # Set Dual Vert Centroid and Edges
        dual_vert.point = Point(x_centroid_avg, y_centroid_avg)
        dual_vert.cc_edges = loop_edges # Dual edges do not necessarily have their other end at thsi point, so they cannot be sorted until the dual is built
        self.dualVerts.add(dual_vert)

    def generate_spatial_lookup(self, resolution: float) -> None:
        self.spatial_lookup_res = resolution
        self.spatial_lookup_ref = []
        for i in range((self.upperXBound - self.lowerXBound) // self.spatial_lookup_res + 1):
            self.spatial_lookup_ref.append([])
            for j in range((self.upperYBound - self.lowerYBound) // self.spatial_lookup_res + 1):
                self.spatial_lookup_ref[j].append([])

        for vert in self.primalVerts:
            x_idx = vert.point.x // resolution
            y_idx = vert.point.y // resolution

            self.spatial_lookup_ref[x_idx][y_idx].append(vert)

    def get_verts_within_radius(self, src: Point, radius: float, role: TwinGraph.VertRole) -> [TwinGraph.Vert]:
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

    def get_closest_vert(self, src: Point, role: TwinGraph.VertRole) -> [TwinGraph.Vert]:
        for radius in range(max(self.upperXBound - self.lowerXBound, self.upperYBound - self.lowerYBound)):
            close = self.get_verts_within_radius(src, radius, role)
            if len(close) > 0:
                closest = min(close, key=lambda vert: Point.dist(src, vert.point))
                return closest

        return None

    class EdgeDir(Enum):
        AB = 1
        BA = 2

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
        cc_edges: [TwinGraph.QuadEdge]

        # Counterclockwise edges must be sorted
        def __init__(self, point: point, weight: float, role: TwinGraph.VertRole, counterclockwiseEdges: [TwinGraph.QuadEdge]=[]) -> None:
            # TODO: Add assert for edge sorting
            self.point = point
            self.weight = weight
            self.role = role
            self.cc_edges = counterclockwiseEdges
            self.link_edges()

        # Edges must already have links to calling vertex
        def register_edges(self, edges: [TwinGraph.QuadEdge]) -> None:
            # TODO: Add assert for edge linking
            edges = self.cc_edges + edges # TODO: Sorted insertion would be faster than re-sorting all
            if self.role.is_primal():
                edges.sort(key=lambda edge: edge.get_primal_rad_from(self)[0])
            if self.role.is_dual():
                edges.sort(key=lambda edge: edge.get_dual_rad_from(self)[0])

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
        primal_annotation: int
        primal_A_cc_next: TwinGraph.QuadEdge
        primal_A_cc_prev: TwinGraph.QuadEdge
        primal_B_cc_next: TwinGraph.QuadEdge
        primal_B_cc_prev: TwinGraph.QuadEdge

        dual_AB: TwinGraph.Vert
        dual_BA: TwinGraph.Vert
        dual_annotation: int
        dual_AB_cc_next: TwinGraph.QuadEdge
        dual_AB_cc_prev: TwinGraph.QuadEdge
        dual_BA_cc_next: TwinGraph.QuadEdge
        dual_BA_cc_prev: TwinGraph.QuadEdge

        def __init__(
                self,
                primal_A: TwinGraph.Vert,
                primal_B: TwinGraph.Vert
            ) -> None:
            self.primal_A = primal_A
            self.primal_B = primal_B
            self.primal_annotation = None
            self.dual_AB = None
            self.dual_BA = None
            self.dual_annotation = None

        # Get the ordered ends of an edge based on an edge direction
        def get_primal_vert_pair(self, dir: TwinGraph.EdgeDir) -> (TwinGraph.Vert, TwinGraph.Vert):
            if dir == TwinGraph.EdgeDir.AB:
                return (self.primal_A, self.primal_B)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.primal_B, self.primal_A)
        def get_dual_vert_pair(self, dir: TwinGraph.EdgeDir) -> (TwinGraph.Vert, TwinGraph.Vert):
            if dir == TwinGraph.EdgeDir.AB:
                return (self.dual_AB, self.dual_BA)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.dual_BA, self.dual_AB)

        # Get the vertex on the other end of an edge along with directional annotation
        def get_primal_dest_from(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, TwinGraph.EdgeDir):
            if src == self.primal_A:
                return (self.primal_B, TwinGraph.EdgeDir.AB)
            elif src == self.primal_B:
                return (self.primal_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_primal_dest_from is not on edge.")
        def get_dual_dest_from(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, TwinGraph.EdgeDir):
            if src == self.dual_AB:
                return (self.dual_BA, TwinGraph.EdgeDir.AB)
            elif src == self.dual_BA:
                return (self.dual_AB, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_dual_dest_from is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def get_primal_rad_from(self, src: TwinGraph.Vert) -> (float, TwinGraph.EdgeDir):
            dest, dir = self.get_primal_dest_from(src)
            return (Point.srcDestRad(src.point, dest.point), dir)
        def get_dual_rad_from(self, src: TwinGraph.Vert) -> (float, TwinGraph.EdgeDir):
            dest, dir = self.get_dual_dest_from(src)
            return (Point.srcDestRad(src.point, dest.point), dir)
        
        # Get the cc next edge from vert
        def get_primal_cc_next_edge(self, vert: TwinGraph.Vert) -> (TwinGraph.QuadEdge, TwinGraph.EdgeDir):
            if vert == self.primal_A:
                _, dir = self.primal_A_cc_next.get_primal_dest_from(vert)
                return (self.primal_A_cc_next, dir)
            if vert == self.primal_B:
                _, dir = self.primal_B_cc_next.get_primal_dest_from(vert)
                return (self.primal_B_cc_next, dir)
        def get_dual_cc_next_edge(self, vert: TwinGraph.Vert) -> (TwinGraph.QuadEdge, TwinGraph.EdgeDir):
            if vert == self.dual_AB:
                _, dir = self.dual_AB_cc_next.get_dual_dest_from(vert)
                return (self.dual_AB_cc_next, dir)
            if vert == self.dual_BA:
                _, dir = self.dual_BA_cc_next.get_dual_dest_from(vert)
                return (self.dual_BA_cc_next, dir)
            
        # Get the cc prev edge from vert
        def get_primal_cc_prev_edge(self, vert: TwinGraph.Vert) -> (TwinGraph.QuadEdge, TwinGraph.EdgeDir):
            if vert == self.primal_A:
                _, dir = self.primal_A_cc_prev.get_primal_dest_from(vert)
                return (self.primal_A_cc_prev, dir)
            if vert == self.primal_B:
                _, dir = self.primal_B_cc_prev.get_primal_dest_from(vert)
                return (self.primal_B_cc_prev, dir)
        def get_dual_cc_prev_edge(self, vert: TwinGraph.Vert) -> (TwinGraph.QuadEdge, TwinGraph.EdgeDir):
            if vert == self.dual_AB:
                _, dir = self.dual_AB_cc_prev.get_dual_dest_from(vert)
                return (self.dual_AB_cc_prev, dir)
            if vert == self.dual_BA:
                _, dir = self.dual_BA_cc_prev.get_dual_dest_from(vert)
                return (self.dual_BA_cc_prev, dir)