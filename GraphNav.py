from __future__ import annotations
from typing import List, Tuple, Set, Optional

import random # TODO: Confirm pseudo-randomness is acceptable

from TwinGraph import *
from RegionTree import *

class GraphNav:
    graph: TwinGraph

    tree_verts: Set[TwinGraph.Vert]
    tree_edges: Set[TwinGraph.QuadEdge]
    bridge_edges: Set[TwinGraph.QuadEdge] # Edges that, if removed, would disconnect the primal tree

    region_tree: RegionTree

    animating: bool
    animation_tracks: List[List[List[Tuple[TwinGraph.QuadEdge, TwinGraph.VertRole, TwinGraph.EdgeDir, int]]]]

    def __init__(self, graph: TwinGraph, region_tree: RegionTree):
        self.graph = graph

        self.tree_verts = set()
        self.tree_edges = set()
        self.bridge_edges = set()

        self.region_tree = region_tree

        self.animating = True
        self.animation_tracks = [[[]], [[]]] # One track for walk, one for region finding

    def walk_division_from(self, region_in: Optional[RegionTree.Region], start_vert: TwinGraph.Vert):
        """
        Perform a walk and region tree update
        Start vert must be withing the region provided
        """
        # Allow inefficient automatic region detection if none provided
        region = None
        if region_in is None:
            for candidate_region in self.region_tree.regions:
                print("Checking region", candidate_region.id_str)
                if start_vert in candidate_region.get_all_dual_verts_contained_by_perimeter(self.graph):
                    region = candidate_region
                    break
            if region is None:
                raise ValueError("Starting vert is not contained in any region.")
        else:
            region = region_in

        # TODO:
         # Verify cc linked list is valid
         # Progressive edge selection logic in traverse_clockwise_loop

        # Build walk than progressively discover regions
        region_edge = self.loop_erased_random_walk_from(start_vert)
        # self.traverse_clockwise_loop(region_edge, TwinGraph.EdgeDir.AB)
        self.develop_region(region, region_edge, TwinGraph.EdgeDir.AB) # Ok to start in arbitrary direction because loop is starting from boundary rather than bridge, so equivalent connections occur from either side
        print("Region Count:", len(self.region_tree.regions))
        self.region_tree.remove_region(region) # Remove old region

        for region in self.region_tree.regions:
            print("Region", region.id_str, "weight", region.weight, "bridges", len(region.bridge_set), "edges", len(region.dual_perimeter))

    def loop_erased_random_walk_from(self, vert: TwinGraph.Vert) -> TwinGraph.QuadEdge:
        assert vert not in self.tree_verts, "Starting vert is already in tree."

        # Check vert is dual
        if vert.role != TwinGraph.VertRole.DUAL:
            raise ValueError("Starting vert role does not match graph_selection.")

        # Walk state
        current_edge: TwinGraph.QuadEdge
        current_vert: TwinGraph.Vert
        consumed_verts: Set[TwinGraph.Vert] = set() # Track for self-collisions
        walk_edges: List[TwinGraph.QuadEdge] = []

        # Set up first edge and vert
        current_edge = vert.cc_edges[random.randint(0, len(vert.cc_edges)-1)]
        current_vert, _ = current_edge.get_dual_dest_from(vert)
        consumed_verts.add(vert)

        while True:
            # Erase loop if present
            if current_vert in consumed_verts:
                rewind_vert: TwinGraph.Vert
                rewind_vert, _ = current_edge.get_dual_dest_from(current_vert)

                while rewind_vert != current_vert: # Since current_vert is in consumed_verts, this will always terminate
                    consumed_verts.remove(rewind_vert)
                    last_edge: TwinGraph.QuadEdge = walk_edges.pop()
                    rewind_vert, _ = last_edge.get_dual_dest_from(rewind_vert)

            # Register progress
            else:
                consumed_verts.add(current_vert)
                walk_edges.append(current_edge)

            # Animation tracking
            if self.animating:
                # Convert walk_edges to the expected tuple format
                walk_tuples = [
                    (edge, TwinGraph.VertRole.DUAL, TwinGraph.EdgeDir.AB, 0)  # Replace FORWARD with actual direction if needed
                    for edge in walk_edges
                ]
                tree_tuples = [
                    (edge, TwinGraph.VertRole.DUAL, TwinGraph.EdgeDir.AB, 0)  # Replace FORWARD with actual direction if needed
                    for edge in self.tree_edges
                ]
                self.animation_tracks[0].append(walk_tuples + tree_tuples)
                
                # Draw bridges
                self.animation_tracks[0][-1].extend([
                    (edge, TwinGraph.VertRole.DUAL, TwinGraph.EdgeDir.AB, 1)
                    for edge in self.bridge_edges
                ])

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
            current_vert, _ = current_edge.get_dual_dest_from(current_vert)
            
        # Commit walk to tree
        for edge in walk_edges:
            self.tree_edges.add(edge)
        for vert in consumed_verts:
            self.tree_verts.add(vert)

        # Detect bridge edges
        for vert in consumed_verts:
            for edge in vert.cc_edges:
                if edge in self.tree_edges:
                    continue
                dest, _ = edge.get_dual_dest_from(vert)
                if dest in self.tree_verts:# and vert.role != TwinGraph.VertRole.DUAL_EXTERIOR and dest.role != TwinGraph.VertRole.DUAL_EXTERIOR:
                    self.bridge_edges.add(edge)
        if self.animating:
            self.animation_tracks[0][-1].extend([
                (edge, TwinGraph.VertRole.DUAL, TwinGraph.EdgeDir.AB, 1)
                for edge in self.bridge_edges
            ])

        # Return edge that began the walk (which is gauranteed to be internal to the region)
        return walk_edges[0]

    def develop_region(self, src_region: RegionTree.Region, start_edge: TwinGraph.QuadEdge, start_dir: TwinGraph.EdgeDir) -> RegionTree.Region:
        # Recursively perform clockwise loops from the start_edge, when a loop crosses a 
        # bridge edge, start a new clockwise loop from the other side of the bridge edge 
        # if the bridge edge is in new doors.  With every completed loop, create a new 
        # region tree vertex that links across the bridge edge to a neighboring region
        # tree vertex).  Build this so as to replace the region passed in.  Build new
        # region tree vertices for each loop found, starting from start_edge

        print("developing region")

        loop = self.traverse_clockwise_loop(start_edge, start_dir)
        root_region = RegionTree.Region(0, loop) # TODO: Get weight of new region
        self.region_tree.add_region(root_region) # Automatically handles region registration

        print("Loop length:", len(loop))

        for edge, dir in loop:
            if edge in self.bridge_edges and edge not in src_region.bridge_set and edge != start_edge:
                # Start a new clockwise loop from the other side of the bridge edge
                neighbor_region = self.develop_region(src_region, edge, dir.reverse())

                # Link the new region tree vertex to the neighbor region tree vertex
                new_region_edge = RegionTree.Edge()
                new_region_edge.end_A = root_region
                new_region_edge.end_B = neighbor_region
                new_region_edge.twin_graph_edge = edge

                self.region_tree.add_edge(new_region_edge) # Automatically handles edge registration

            if edge in self.bridge_edges and edge in src_region.bridge_set and edge != start_edge:
                # Edge is a bridge to src_region, replace old link to src_region with new link to root_region
                old_region_edge = src_region.bridge_to_region_edge_map[edge]
                other_region, _ = old_region_edge.get_dest_from(src_region)
                # Not necessary to unregister old edge because region is being deleted

                new_region_edge = RegionTree.Edge()
                new_region_edge.end_A = root_region
                new_region_edge.end_B = other_region
                new_region_edge.twin_graph_edge = edge

                self.region_tree.add_edge(new_region_edge) # Automatically handles edge registration

        return root_region

    def traverse_clockwise_loop(self, start_edge: TwinGraph.QuadEdge, edge_dir: TwinGraph.EdgeDir) -> List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]:
        """
        Traverses a clockwise loop within the tree and bridge edges, starting from start_edge.
        Returns the list of visited edges in order.
        """

        visited_edges: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] = []

        # Determine starting vert from edge and direction
        current_vert, _ = start_edge.get_dual_vert_pair(edge_dir)
        current_edge = start_edge

        # Mark start vert
        start_vert = current_vert

        print("Starting loop at vert", current_vert.id_str, "on edge", current_edge.id_str, "going", edge_dir.name)
        while True:
            # Move to the next vert along the current edge and direction
            current_vert, current_dir = current_edge.get_dual_dest_from(current_vert)
            visited_edges.append((current_edge, current_dir))

            # if self.animating:
            #     # Convert visited_edges to the expected tuple format
            #     traverse_link = (current_edge, TwinGraph.VertRole.DUAL, current_dir, 2) 
            #     self.animation_tracks[0].append(self.animation_tracks[0][-1] + [traverse_link])

            # Find the next edge ccw that is on the tree or a bridge
            next_edge = current_edge
            next_dir = current_dir
            while True:
                if current_vert == next_edge.dual_BA:
                    next_edge = next_edge.dual_BA_cc_next
                else:
                    next_edge = next_edge.dual_AB_cc_next
                _, next_dir = next_edge.get_dual_dest_from(current_vert)

                print("At vert", current_vert.id_str, "checking edge to", next_edge.get_dual_dest_from(current_vert)[0].id_str, "along", next_edge.get_dual_rad_from(current_vert)[0])

                if next_edge in self.tree_edges or next_edge in self.bridge_edges:
                    break
            current_edge = next_edge
            current_dir = next_dir

            # Check if we've looped back to the start
            if current_edge == start_edge and current_vert == start_vert:
                break

        return visited_edges
            
