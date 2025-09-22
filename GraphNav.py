from __future__ import annotations
from typing import List, Tuple, Set, Optional

import random # TODO: Confirm pseudo-randomness is acceptable

from TwinGraph import *

class GraphNav:
    type GraphSelection = TwinGraph.VertRole

    graph: TwinGraph
    graph_selection: GraphSelection
    tree_verts: Set[TwinGraph.Vert]
    tree_edges: Set[TwinGraph.Vert]

    animating: bool
    animation_track: List[List[Tuple[TwinGraph.QuadEdge, TwinGraph.VertRole, TwinGraph.EdgeDir]]]

    def __init__(self, graph: TwinGraph, graph_selection: GraphSelection):
        self.graph = graph
        self.graph_selection = graph_selection

        self.animating = True
        self.animation_track = []

    def loop_erased_random_walk_from(self, vert: TwinGraph.Vert):
        current_vert: TwinGraph.Vert = vert
        current_edge: TwinGraph.QuadEdge = vert.cc_edges[random.randint(0, len(vert.cc_edges)-1)]
        consumed_verts: Set[TwinGraph.Vert] = set() # Track for self-collisions
        walk_edges: List[TwinGraph.QuadEdge] = []

        while True:
            # Erase loop if present
            if current_vert in consumed_verts:
                rewind_vert: Optional[TwinGraph.Vert] = None
                rewind_vert = current_edge.get_primal_dest_from(current_vert)

                while rewind_vert != current_vert:
                    last_edge: TwinGraph.QuadEdge = walk_edges.pop()
                    rewind_vert = last_edge.get_primal_dest_from(rewind_vert)

            # Register progress

            # Find next vert and edge
            current_edge = current_vert.cc_edges[random.randint(0, len(vert.cc_edges)-1)]
            current_vert, _ = current_edge.get_primal_dest_from(current_vert)
            
