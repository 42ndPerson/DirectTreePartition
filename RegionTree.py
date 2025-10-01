from __future__ import annotations
from typing import List, Tuple, Set

from TwinGraph import *

class RegionTree:
    verts: Set[RegionTree.Vert]
    edges: Set[RegionTree.Edge]

    type EdgeDir = TwinGraph.EdgeDir

    def __init__(self, total_weight: float) -> None:
        self.verts = {RegionTree.Vert(total_weight)}
        self.edges = Set()

    class Vert:
        weight: float
        cc_edges: List[RegionTree.Edge]

        def __init__(self, weight: float, edges: List[RegionTree.Edge]=[]) -> None:
            self.weight = weight

        def register_edges(self, edges: List[RegionTree.Edge]) -> None:
            # TODO: Add assert for edge linking
            edges = self.cc_edges + edges # TODO: Sorted insertion would be faster than re-sorting all
            edges.sort(key=lambda edge: edge.get_rad_from(self)[0])

            self.cc_edges = edges
            # TODO: Add edge linking if needed

    class Edge:
        # Edge A and B ends must match directionally with underlying twin_graph_edge
        end_A: RegionTree.Vert
        end_B: RegionTree.Vert
        twin_graph_edge: TwinGraph.QuadEdge

        # Get the vertex on the other end of an edge along with directional annotation
        def get_dest_from(self, src: TwinGraph.Vert) -> Tuple[TwinGraph.Vert, TwinGraph.EdgeDir]:
            if src == self.end_A:
                return (self.end_B, TwinGraph.EdgeDir.AB)
            elif src == self.end_B:
                return (self.end_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_primal_dest_from is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation  
        def get_rad_from(self, src: TwinGraph.Vert) -> Tuple[float, TwinGraph.EdgeDir]:
            _, dir = self.get_dest_from(src)
            return (self.twin_graph_edge.get_primal_rad_along(dir), dir)
