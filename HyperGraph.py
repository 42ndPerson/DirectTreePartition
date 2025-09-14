from __future__ import annotations
from enum import Enum
from Euclid import *

class HyperGraph:
    primalVerts: set(HyperGraph.HyperVert)
    dualVerts: set(HyperGraph.HyperVert)
    edges: set(HyperGraph.HyperEdge)

    lowerXBound: float
    upperXBound: float
    lowerYBound: float
    upperYBound: float

    spatial_lookup_res: float
    spatial_lookup_ref: list[list[HyperGraph.HyperVert]] # Indexing on x and then y

    def __init__(self, points: [Point], weights: [float], edgeIdxs: [(int,int)]) -> None:
        assert len(points) == len(weights)
        zippedPointWeights = zip(points, weights)

        # Create key data structures
        self.primalVerts = [HyperGraph.HyperVert(point, weight, HyperGraph.VertRole.PRIMAL) for point, weight in zippedPointWeights]
        self.edges = [HyperGraph.HyperEdge(self.primalVerts[idxA], self.primalVerts[idxB]) for idxA, idxB in edgeIdxs]

        # Note bounds
        self.lowerXBound = min(points, key=lambda a: a.x).x
        self.upperXBound = max(points, key=lambda a: a.x).x
        self.lowerYBound = min(points, key=lambda a: a.y).y
        self.upperYBound = max(points, key=lambda a: a.y).y

        # Register edges with vertices
        vertEdgeDict = {}
        for edge in self.edges:
            vertA = edge.primalA
            vertB = edge.primalB
            if vertA not in vertEdgeDict:
                vertEdgeDict[vertA] = []
            if vertB not in vertEdgeDict:
                vertEdgeDict[vertB] = []
            vertEdgeDict[vertA].append(edge)
            vertEdgeDict[vertB].append(edge)
        for vert, vertEdges in vertEdgeDict.items():
            vert.addEdges(vertEdges)

        # Construct dual graph
        # faceRoot = self.primalVerts[0]
            
    def construct_dual(self) -> None:
        consumed_edges_AB: set(HyperGraph.HyperEdge) = set()
        consumed_edges_BA: set(HyperGraph.HyperEdge) = set()
        
        face_root = self.primalVerts[0]
        total_angle = 0

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

    def get_verts_within_radius(self, src: Point, radius: float, role: HyperGraph.VertRole) -> [HyperGraph.HyperVert]:
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

    def get_closest_vert(self, src: Point, role: HyperGraph.VertRole) -> [HyperGraph.HyperVert]:
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

    class EdgeDir(Enum):
        AB = 1
        BA = -1

        def op(self) -> HyperGraph.EdgeDir:
            return HyperGraph.EdgeDir(-self.value)

    class HyperVert:
        point: Point
        weight: float
        role: HyperGraph.VertRole
        counterclockwiseEdges: [HyperGraph.HyperEdge]

        def __init__(self, point: point, weight: float, role: HyperGraph.VertRole) -> None:
            self.point = point
            self.weight = weight
            self.role = role
            self.clockwiseEdges = []

        def addEdges(self, edges: [HyperGraph.HyperEdge]) -> None:
            edges = self.clockwiseEdges + edges
            if self.role.is_primal():
                edges.sort(key=lambda edge: edge.getPrimalRadFrom(self)[0])
            if self.role.is_dual():
                edges.sort(key=lambda edge: edge.getDualRadFrom(self)[0])

            self.counterclockwiseEdges = edges

    # HyperEdge contains both primal and dual edge info
    class HyperEdge:
        primalA: HyperGraph.HyperVert
        primalB: HyperGraph.HyperVert
        primalAnnotation: int
        dualA: HyperGraph.HyperVert
        dualB: HyperGraph.HyperVert
        dualAnnotation: int

        def __init__(
                self,
                primalA: HyperGraph.HyperVert,
                primalB: HyperGraph.HyperVert
            ) -> None:
            self.primalA = primalA
            self.primalB = primalB
            self.primalAnnotation = None
            self.dualAB = None
            self.dualBA = None
            self.dualAnnotation = None

        # Get the vertex on the other end of an edge along with directional annotation
        def getPrimalDestFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, HyperGraph.EdgeDir):
            if src == self.primalA:
                return (self.primalB, HyperGraph.EdgeDir.AB)
            elif src == self.primalB:
                return (self.primalA, HyperGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for getPrimalDestFrom is not on edge.")
        def getDualDestFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, HyperGraph.EdgeDir):
            if src == self.dualA:
                return (self.dualB, HyperGraph.EdgeDir.AB)
            elif src == self.dualB:
                return (self.dualA, HyperGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for getDualDestFrom is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def getPrimalRadFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, float):
            dest, dir = self.getPrimalDestFrom(src)
            return (Point.srcDestRad(src.point, dest.point), dir)
        def getDualRadFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, float):
            dest, dir = self.getDualDestFrom(src)
            return (Point.srcDestRad(src.point, dest.point), dir)