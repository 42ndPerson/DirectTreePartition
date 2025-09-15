from __future__ import annotations

from collections import deque
import bisect

from enum import Enum
from Euclid import *

class TwinGraph:
    primalVerts: set(TwinGraph.Vert)
    dualVerts: set(TwinGraph.Vert)
    edges: set(TwinGraph.HyperEdge)

    lowerXBound: float
    upperXBound: float
    lowerYBound: float
    upperYBound: float

    spatial_lookup_res: float
    spatial_lookup_ref: [[TwinGraph.Vert]] # Indexing on x and then y

    animating: bool
    animation_frame: int
    animation_track: [[(TwinGraph.HyperEdge, TwinGraph.VertRole, TwinGraph.EdgeDir)]]

    def __init__(self, points: [Point], weights: [float], edgeIdxs: [(int,int)]) -> None:
        assert len(points) == len(weights)
        zippedPointWeights = zip(points, weights)

        # Create key data structures
        self.primalVerts = [TwinGraph.Vert(point, weight, TwinGraph.VertRole.PRIMAL) for point, weight in zippedPointWeights]
        self.edges = [TwinGraph.HyperEdge(self.primalVerts[idxA], self.primalVerts[idxB]) for idxA, idxB in edgeIdxs]

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
            vert.addEdges(vertEdges)

        # Aux Data
        # Animation
        self.animating = True
        self.animation_frame = 0
        self.animation_track = []

        # Construct dual graph
        # faceRoot = self.primalVerts[0]
            
    def construct_dual(self) -> None:
        vert_next_outgoing_edge: {TwinGraph.Vert: int} = {}
        edge_deque: deque((TwinGraph.HyperEdge, TwinGraph.EdgeDir)) = deque()
        consumed_edges_AB: set(TwinGraph.HyperEdge) = set()
        consumed_edges_BA: set(TwinGraph.HyperEdge) = set()

        edge_deque.append(next(iter(self.edges, TwinGraph.EdgeDir.AB)))

        while len(edge_deque) > 0:
            edge, _ = edge_deque.popleft()
            self.construct_face_loop(
                edge, 
                dir, 
                edge_deque, 
                vert_next_outgoing_edge, 
                consumed_edges_AB, 
                consumed_edges_BA
            )

    def construct_face_loop(
            self,
            src_edge: TwinGraph.HyperEdge,
            src_dir: TwinGraph.EdgeDir,
            edge_deque: deque((TwinGraph.HyperEdge, TwinGraph.EdgeDir)),
            vert_next_outgoing_edge: {TwinGraph.Vert: int},
            consumed_edges_AB: set(TwinGraph.HyperEdge),
            consumed_edges_BA: set(TwinGraph.HyperEdge)):
        
        edges: [TwinGraph.HyperEdge] = []
        
        origin_vert, _ = src_edge.get_primal_vert_pair(src_dir)

        current_edge: TwinGraph.HyperEdge = src_edge
        while True:
            pass

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
        cc_edges: [TwinGraph.HyperEdge]

        def __init__(self, point: point, weight: float, role: TwinGraph.VertRole) -> None:
            self.point = point
            self.weight = weight
            self.role = role
            self.cc_edges = []

        # Counterclockwise edges must be sorted
        def __init__(self, point: point, weight: float, role: TwinGraph.VertRole, counterclockwiseEdges: [TwinGraph.HyperEdge]) -> None:
            self.point = point
            self.weight = weight
            self.role = role
            self.cc_edges = counterclockwiseEdges

        def addEdges(self, edges: [TwinGraph.HyperEdge]) -> None:
            edges = self.clockwiseEdges + edges
            if self.role.is_primal():
                edges.sort(key=lambda edge: edge.getPrimalRadFrom(self)[0])
            if self.role.is_dual():
                edges.sort(key=lambda edge: edge.getDualRadFrom(self)[0])

            for i, edge in enumerate(edges):
                if self.role.is_primal():
                    edges.sort(key=lambda edge: edge.getPrimalRadFrom(self)[0])
                if self.role.is_dual():
                    edges.sort(key=lambda edge: edge.getDualRadFrom(self)[0])
            self.cc_edges = edges

    class EdgePair:
        edge: TwinGraph.DirectedEdge

    class EdgeNexus:
        primal_A: TwinGraph.DirectedEdge
        primal_B: TwinGraph.DirectedEdge
        dual_A: TwinGraph.DirectedEdge
        dual_B: TwinGraph.DirectedEdge

    class DirectedEdge:
        src: TwinGraph.Vert
        dest: TwinGraph.Vert
        nexus: TwinGraph.EdgeNexus
        annotation: int

        src_cc_prev: TwinGraph.DirectedEdge
        src_cc_next: TwinGraph.DirectedEdge
        prev_cc_prev: TwinGraph.DirectedEdge
        prev_cc_next: TwinGraph.DirectedEdge

        def __init__(self, src: TwinGraph.Vert, dest: TwinGraph.Vert):
            self.src = src
            self.dest = dest

            self.annotation = None
            self.src_cc_prev = None
            self.src_cc_next = None
            self.dest_cc_prev = None
            self.dest_cc_next = None

    # HyperEdge contains both primal and dual edge info
    class HyperEdge:
        primal_A: TwinGraph.Vert
        primal_B: TwinGraph.Vert
        primal_annotation: int
        primal_A_cc_next: TwinGraph.HyperEdge
        primal_A_cc_prev: TwinGraph.HyperEdge
        primal_B_cc_next: TwinGraph.HyperEdge
        primal_B_cc_prev: TwinGraph.HyperEdge

        dual_A: TwinGraph.Vert
        dual_B: TwinGraph.Vert
        dual_annotation: int
        dual_A_cc_next: TwinGraph.HyperEdge
        dual_A_cc_prev: TwinGraph.HyperEdge
        dual_B_cc_next: TwinGraph.HyperEdge
        dual_B_cc_prev: TwinGraph.HyperEdge

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
                return (self.dual_A, self.dual_B)
            if dir == TwinGraph.EdgeDir.BA:
                return (self.dual_B, self.dual_A)

        # Get the vertex on the other end of an edge along with directional annotation
        def getPrimalDestFrom(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, TwinGraph.EdgeDir):
            if src == self.primal_A:
                return (self.primal_B, TwinGraph.EdgeDir.AB)
            elif src == self.primal_B:
                return (self.primal_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for getPrimalDestFrom is not on edge.")
        def getDualDestFrom(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, TwinGraph.EdgeDir):
            if src == self.dual_A:
                return (self.dual_B, TwinGraph.EdgeDir.AB)
            elif src == self.dual_B:
                return (self.dual_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for getDualDestFrom is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def getPrimalRadFrom(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, float):
            dest, dir = self.getPrimalDestFrom(src)
            return (Point.srcDestRad(src.point, dest.point), dir)
        def getDualRadFrom(self, src: TwinGraph.Vert) -> (TwinGraph.Vert, float):
            dest, dir = self.getDualDestFrom(src)
            return (Point.srcDestRad(src.point, dest.point), dir)