from __future__ import annotations
from typing import List, Tuple, Set, Optional

import random # TODO: Confirm pseudo-randomness is acceptable

from TwinGraph import *

class GraphNav:
    type GraphSelection = TwinGraph.VertRole

    graph: TwinGraph
    graph_selection: GraphSelection
    tree_verts: Set[TwinGraph.Vert]
    tree_edges: Set[TwinGraph.QuadEdge]

    animating: bool
    animation_tracks: List[List[List[Tuple[TwinGraph.QuadEdge, TwinGraph.VertRole, TwinGraph.EdgeDir]]]]

    def __init__(self, graph: TwinGraph, graph_selection: GraphSelection):
        self.graph = graph
        self.graph_selection = graph_selection

        self.tree_verts = set()
        self.tree_edges = set()

        self.animating = True
        self.animation_tracks = [[[]]]

    def loop_erased_random_walk_from(self, vert: TwinGraph.Vert):
        print("Walking")
        # Check vert is of role graph_selection
        if vert.role != self.graph_selection:
            raise ValueError("Starting vert role does not match graph_selection.")

        # Walk state
        current_edge: TwinGraph.QuadEdge
        current_vert: TwinGraph.Vert
        consumed_verts: Set[TwinGraph.Vert] = set() # Track for self-collisions
        walk_edges: List[TwinGraph.QuadEdge] = []

        # Set up first edge and vert
        current_edge = vert.cc_edges[random.randint(0, len(vert.cc_edges)-1)]
        current_vert, _ = (
            current_edge.get_primal_dest_from(vert) if 
            self.graph_selection == TwinGraph.VertRole.PRIMAL else 
            current_edge.get_dual_dest_from(vert)
        )
        consumed_verts.add(vert)

        while True:
            # Erase loop if present
            if current_vert in consumed_verts:
                rewind_vert: TwinGraph.Vert
                rewind_vert, _ = (
                    current_edge.get_primal_dest_from(current_vert) if 
                    self.graph_selection == TwinGraph.VertRole.PRIMAL else 
                    current_edge.get_dual_dest_from(current_vert)
                )

                while rewind_vert != current_vert: # Since current_vert is in consumed_verts, this will always terminate
                    consumed_verts.remove(rewind_vert)
                    last_edge: TwinGraph.QuadEdge = walk_edges.pop()
                    rewind_vert, _ = (
                        last_edge.get_primal_dest_from(rewind_vert) if 
                        self.graph_selection == TwinGraph.VertRole.PRIMAL else 
                        last_edge.get_dual_dest_from(rewind_vert)
                    )

            # Register progress
            else:
                consumed_verts.add(current_vert)
                walk_edges.append(current_edge)

            # Animation tracking
            if self.animating:
                # Convert walk_edges to the expected tuple format
                walk_tuples = [
                    (edge, self.graph_selection, TwinGraph.EdgeDir.AB)  # Replace FORWARD with actual direction if needed
                    for edge in walk_edges
                ]
                tree_tuples = [
                    (edge, self.graph_selection, TwinGraph.EdgeDir.AB)  # Replace FORWARD with actual direction if needed
                    for edge in self.tree_edges
                ]
                self.animation_tracks[0].append(walk_tuples + tree_tuples)

            # Check if done
            if (
                current_vert in self.tree_verts or # Walk hits tree
                current_vert.role == TwinGraph.VertRole.DUAL_EXTERIOR or # Dual walk hits exterior
                ((current_edge.dual_AB is not None and current_edge.dual_AB.role == TwinGraph.VertRole.DUAL_EXTERIOR) or # Primal walk moves along edge with exterior dual
                 (current_edge.dual_BA is not None and current_edge.dual_BA.role == TwinGraph.VertRole.DUAL_EXTERIOR))
            ):
                break

            # Find next vert and edge
            current_edge = current_vert.cc_edges[random.randint(0, len(current_vert.cc_edges)-1)]
            current_vert, _ = (
                    current_edge.get_primal_dest_from(current_vert) if 
                    self.graph_selection == TwinGraph.VertRole.PRIMAL else 
                    current_edge.get_dual_dest_from(current_vert)
                )
            
        for edge in walk_edges:
            self.tree_edges.add(edge)
        for vert in consumed_verts:
            self.tree_verts.add(vert)
            
