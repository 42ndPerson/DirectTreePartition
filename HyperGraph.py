from enum import Enum
import Euclid

class HyperGraph:
    primalVerts: set(HyperGraph.HyperVert)
    dualVerts: set(HyperGraph.HyperVert)
    edges: set(HyperGraph.HyperEdge)

    def __init__(self, points: [Euclid.Point], weights: [float], edges: [(int,int)]) -> None:
        assert points.length == weights.length
        zippedPointWeights = zip(points, weights)

        # Create key data structures
        self.primalVerts = map(lambda point, weight: HyperGraph.HyperVert(point, weight), zippedPointWeights)
        self.edges = map(lambda idxA, idxB: HyperGraph.HyperEdge(self.primalVerts[idxA], self.primalVerts[idxB]), edges)

        # Register edges with verticies
        vertEdgeDict = {}
        for edge in self.edges:
            vert = edge.primalA
            if vert not in vertEdgeDict:
                vertEdgeDict[vert] = []
            vertEdgeDict[vert].append(edge)
        for vert, edges in vertEdgeDict.items():
            vert.addEdges(edges)

        # Construct dual graph
        # faceRoot = self.primalVerts[0]



    class VertRole(Enum):
        PRIMAL = 1
        DUAL = 2

    class EdgeDir(Enum):
        AB = 1
        BA = -1

        def op(self) -> HyperGraph.EdgeDir:
            return HyperGraph.EdgeDir(-self.value)

    class HyperVert:
        point: Euclid.Point
        weight: float
        role: HyperGraph.VertRole
        counterclockwiseEdges: [HyperGraph.HyperEdge]

        def __init__(self, point: Euclid.point, weight: float, role: HyperGraph.VertRole) -> None:
            self.point = point
            self.weight = weight
            self.role = role
            self.clockwiseEdges = []

        def addEdges(self, edges: [HyperGraph.HyperEdge]) -> None:
            edges = self.clockwiseEdges + edges
            if self.role == HyperGraph.VertRole.PRIMAL:
                edges.sort(key=lambda edge: edge.getPrimalDestFromRad(self)[0])
            if self.role == HyperGraph.VertRole.DUAL:
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
        def getPrimalDestFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, EdgeDir):
            if src == self.primalA:
                return (self.primalB, EdgeDir.AB)
            elif src == self.primalB:
                return (self.primalA, EdgeDir.BA)
            else:
                raise KeyError("Src for getPrimalDestFrom is not on edge.")
        def getDualDestFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, EdgeDir):
            if src == self.dualA:
                return (self.dualB, EdgeDir.AB)
            elif src == self.dualB:
                return (self.dualA, EdgeDir.BA)
            else:
                raise KeyError("Src for getDualDestFrom is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def getPrimalRadFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, float):
            dest, dir = self.getPrimalDestFrom(src)
            return (Euclid.srcDestRad(src.point, dest.point), dir)
        def getDualRadFrom(self, src: HyperGraph.HyperVert) -> (HyperGraph.HyperVert, float):
            dest, dir = self.getDualDestFrom(src)
            return (Euclid.srcDestRad(src.point, dest.point), dir)