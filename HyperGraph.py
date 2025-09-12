from __future__ import annotations
from enum import Enum
import Euclid

class HyperGraph:
    primalVerts: set(HyperGraph.HyperVert)
    dualVerts: set(HyperGraph.HyperVert)
    edges: set(HyperGraph.HyperEdge)

    class VertType(Enum):
        PRIMAL = 1
        DUAL = 2

    class EdgeDir(Enum):
        AB = 1
        BA = -1

        def op(self):
            return HyperGraph.EdgeDir(-self.value)

    class HyperVert:
        point: Euclid.Point
        weight: float
        role: HyperGraph.VertType
        clockwiseEdges: list[HyperGraph.HyperEdge]

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