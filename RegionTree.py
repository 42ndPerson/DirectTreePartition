from __future__ import annotations
import bisect
from itertools import pairwise
import random
from typing import List, Tuple, Set, NamedTuple
from collections import deque
from enum import Enum
import warnings
from time import sleep
from operator import attrgetter
import numpy as np

from TwinGraph import *
from Euclid import *
import time

# Lightweight tuple for boundary entries used in get_uniform_interior_vert
class BoundaryElementNT(NamedTuple):
    # boundary_type: 0 = PERIMETER, 1 = LINING
    dfs_index: int
    boundary_type: int
    vert: TwinGraph.Vert

class RegionTree:
    __slots__ = (
        'regions', 
        'edges', 
        'graph', 
        'total_weight',
        'epsilon', 
        'edge_centers',  
        'id_str'
    )

    regions: Set[RegionTree.Region]
    edges: Set[RegionTree.Edge]

    graph: TwinGraph

    # Central region and edge
    total_weight: float
    epsilon: float
    edge_centers: Set[RegionTree.Edge]

    # Instance labeling
    instance_counter: int = 0
    id_str: str

    type EdgeDir = TwinGraph.EdgeDir

    def __init__(self, graph: TwinGraph, epsilon: float = 0, src_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]=[]) -> None:
        # The src region of a region tree is perimter-less
        # The only perimeter it could have would be the all
        #  the edges to the dual exterior listed twice. 
        # However, this would immediately cut off any potential
        #  primal trees that connect along the edges of the 
        #  graph.
        # We are able to omit the perimeter to solve this
        #  problem because the graph is finite, so all paths
        #  are contained naturally by the graph structure.
        self.total_weight = graph.total_primal_weight
        self.epsilon = epsilon # Allowable deviation from half of total weight for each half (difference is 2 * epsilon * total_weight)
        
        self.regions = set()
        self.add_region(RegionTree.Region(self.total_weight, src_perimeter)) # Start with single region containing whole graph
        self.edges = set()

        self.graph = graph

        self.edge_centers = set()

        self.id_str = f"RTree{RegionTree.instance_counter}"
        RegionTree.instance_counter += 1

    def add_region(self, vert: RegionTree.Region) -> None:
        self.regions.add(vert) 
        vert.region_tree = self

    def remove_region(self, vert: RegionTree.Region) -> None:
        for edge in list(vert.region_edges):
            self.remove_edge(edge)
        self.regions.remove(vert)

    def add_edge(self, edge: RegionTree.Edge) -> None:
        self.edges.add(edge)
        edge.end_A.register_edge(edge)
        edge.end_B.register_edge(edge)
        edge.region_tree = self

    def remove_edge(self, edge: RegionTree.Edge) -> None:
        self.edges.remove(edge)
        edge.end_A.unregister_edge(edge)
        edge.end_B.unregister_edge(edge)

    def debug_inventory(self, edge: RegionTree.Edge) -> None:
        print("Inventory on Edge:", edge.id_str)
        print("    ", "Edge in region tree" if edge in self.edges else "Edge NOT in region tree")
        print("    ", "Edge end_A in region tree" if edge.end_A in self.regions else "Edge end_A NOT in region tree")
        print("    ", "Edge end_B in region tree" if edge.end_B in self.regions else "Edge end_B NOT in region tree")
        print("    ", "Edge in end_A region edges" if edge in edge.end_A.region_edges else "Edge NOT in end_A region edges")
        print("    ", "Edge in end_B region edges" if edge in edge.end_B.region_edges else "Edge NOT in end_B region edges")
        print("    ", "Edge in end_A bridge_to_region_edge_map" if edge in edge.end_A.bridge_to_region_edge_map.values() else "Edge NOT in end_A bridge_to_region_edge_map")
        print("    ", "Edge in end_B bridge_to_region_edge_map" if edge in edge.end_B.bridge_to_region_edge_map.values() else "Edge NOT in end_B bridge_to_region_edge_map")
        print("    ", "Edge twin_graph_edge in end_A bridge_set" if edge.twin_graph_edge in edge.end_A.bridge_set else "Edge twin_graph_edge NOT in end_A bridge_set")
        print("    ", "Edge twin_graph_edge in end_B bridge_set" if edge.twin_graph_edge in edge.end_B.bridge_set else "Edge twin_graph_edge NOT in end_B bridge_set")
        print("    ", "Edge twin_graph_edge in end_A bridge_to_region_edge_map" if edge.twin_graph_edge in edge.end_A.bridge_to_region_edge_map.keys() else "Edge twin_graph_edge NOT in end_A bridge_to_region_edge_map")
        print("    ", "Edge twin_graph_edge in end_B bridge_to_region_edge_map" if edge.twin_graph_edge in edge.end_B.bridge_to_region_edge_map.keys() else "Edge twin_graph_edge NOT in end_B bridge_to_region_edge_map")

        print("    ", "--------- Scanning all regions for edge presence ---------")
        for region in self.regions:
            if edge in region.region_edges:
                print("    ", "Edge found in region edges of", region.id_str)
            if edge in region.bridge_to_region_edge_map.values():
                print("    ", "Edge found in bridge_to_region_edge_map of", region.id_str)

    class Region:
        __slots__ = (
            'weight', 
            'dual_perimeter',
            'dual_perimeter_edges', 
            'region_edges', 
            'cc_region_edges', 
            'bridge_set', 
            'bridge_to_region_edge_map', 
            'region_tree', 
            'point', 
            'id_str'
        )

        # Attributes
        weight: int # Total weight of the region
        dual_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] # Edges that define the perimeter of the region in the dual graph
        dual_perimeter_edges: Set[TwinGraph.QuadEdge] # Edges that define the perimeter of the region in the dual graph

        # Connections
        region_edges: Set[RegionTree.Edge] # Region tree edges connected to this vertex
        cc_region_edges: List[RegionTree.Edge] # Counter-clockwise ordered region edges around this vertex
        bridge_set: Set[TwinGraph.QuadEdge] # Edges that, if removed, would disconnect the primal tree
        bridge_to_region_edge_map: Dict[TwinGraph.QuadEdge, RegionTree.Edge] # Map from bridge edges to region tree edges

        # Tree
        region_tree: RegionTree # Must be set when region is added to tree

        # Visualization
        point: Point

        # Instance labeling
        instance_counter: int = 0
        id_str: str

        def __init__(self, weight: int, dual_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]) -> None:
            self.weight = weight
            self.dual_perimeter = dual_perimeter
            self.dual_perimeter_edges = {edge for edge, _ in dual_perimeter}

            self.region_edges = set()
            self.bridge_set = set()
            self.bridge_to_region_edge_map = dict()

            self.id_str = f'R{RegionTree.Region.instance_counter:x}' # Hexadecimal ID for easier reading in debug
            RegionTree.Region.instance_counter += 1

            self.calc_point()

        def register_edge(self, edge: RegionTree.Edge) -> None:
            self.region_edges.add(edge)
            self.bridge_set.add(edge.twin_graph_edge)
            self.bridge_to_region_edge_map[edge.twin_graph_edge] = edge

        def unregister_edge(self, edge: RegionTree.Edge) -> None:
            self.region_edges.remove(edge)
            self.bridge_set.remove(edge.twin_graph_edge)
            del self.bridge_to_region_edge_map[edge.twin_graph_edge]

        def generate_cc_region_edges(self) -> None: # This should be resiliant without valid dual vert point placement
            self.cc_region_edges = []
            for edge, _ in reversed(self.dual_perimeter):
                if edge in self.bridge_to_region_edge_map:
                    self.cc_region_edges.append(self.bridge_to_region_edge_map[edge])

        def calc_point(self) -> None:
            x_sum = 0.0
            y_sum = 0.0
            scaling = 0.0
            for edge, dir in self.dual_perimeter:
                src, _ = edge.get_dual_vert_pair(dir)
                if src.role != TwinGraph.VertRole.DUAL_EXTERIOR:
                    x_sum += src.point.x
                    y_sum += src.point.y
                    scaling += 1

            if len(self.dual_perimeter) == 0:
                self.point = Point(0.0, 0.0)
            else:
                self.point = Point(x_sum / scaling, y_sum / scaling)
            
        def get_perimeter_verts(self) -> Set[TwinGraph.Vert]: 
            # Worthwhile to compute on demand as is only called when finding 
            # starting point, which is an operation limited to few regions
            verts = set()
            for edge, dir in self.dual_perimeter:
                src, _ = edge.get_dual_vert_pair(dir)
                verts.add(src)
            return verts

        def get_interior_lining_verts(self) -> Tuple[Set[TwinGraph.Vert], Set[TwinGraph.Vert]]:
            """
            Returns a set of dual verts that line the inside of the perimeter, and
            a set of dual verts that comprise the perimeter.
            """

            if len(self.dual_perimeter) == 0:
                return ({e.get_dual_dest_from(self.region_tree.graph.external_dual_vert)[0] 
                        for e in self.region_tree.graph.external_dual_vert.cc_edges}, 
                        set([self.region_tree.graph.external_dual_vert]))

            full_perimeter_edges = self.dual_perimeter_edges
            lining_verts: Set[TwinGraph.Vert] = set()
            perimeter_verts = set()

            for edge, dir in self.dual_perimeter:
                rotary_center: TwinGraph.Vert
                _, rotary_center = edge.get_dual_vert_pair(dir)
                perimeter_verts.add(rotary_center)

                working_edge = edge
                while True:
                    working_edge, _ = working_edge.get_dual_cc_next_edge(rotary_center)
                    if working_edge in full_perimeter_edges: 
                        break

                    lining_verts.add(working_edge.get_dual_dest_from(rotary_center)[0])
                    # Performance note: Set is absorbing ~50% double vert detections with grid graphs

            return lining_verts, perimeter_verts

        # Performance notes:
#         ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#     73786    0.464    0.000  310.981    0.004 RegionTree.py:233(get_uniform_interior_vert)
#     71596   46.859    0.001  146.198    0.002 RegionTree.py:242(_build_perimeter_table)
#     71596   42.450    0.001  107.536    0.002 RegionTree.py:202(get_interior_lining_verts)
#     71596   28.467    0.000   34.416    0.000 RegionTree.py:269(_collect_dual_tree_path_pairs)
#     71596    0.031    0.000   30.261    0.000 RegionTree.py:266(_perimeter_table_sort)
#   1050494    3.060    0.000   29.630    0.000 RegionTree.py:147(__init__)
#   1050494   17.936    0.000   26.570    0.000 RegionTree.py:177(calc_point)
#     71496   13.922    0.000   15.467    0.000 RegionTree.py:284(_consolidate_ranges_and_count_weight)
#   4277944    2.496    0.000    8.712    0.000 RegionTree.py:80(add_edge)
        # Interior lining verts and the non-sorting parts of build_perimeter_table are the bottlenecks
        def get_uniform_interior_vert(self) -> Optional[TwinGraph.Vert]:
            if len(self.dual_perimeter) == 0:
                return random.choice(list(self.region_tree.graph.dualVerts))

            # Use simple integer tags instead of Enum for performance: 0=PERIMETER, 1=LINING
            PERIMETER = 0
            LINING = 1

            # -------------------- Helper sub-functions --------------------
            def _build_perimeter_table(perimeter_verts: Set[TwinGraph.Vert],
                                       lining_verts: Set[TwinGraph.Vert]) -> List[BoundaryElementNT]:
                perimeter_table: List[BoundaryElementNT] = \
                    [BoundaryElementNT(v.map_tree_vert.dfs_in, PERIMETER, v) for v in perimeter_verts] \
                    + [BoundaryElementNT(v.map_tree_vert.dfs_out, PERIMETER, v) for v in perimeter_verts] \
                    + [BoundaryElementNT(v.map_tree_vert.dfs_in, LINING, v) for v in lining_verts] \
                    + [BoundaryElementNT(v.map_tree_vert.dfs_out, LINING, v) for v in lining_verts]
                # Sort boundary elements by dfs index to group path pairs
                # Use a key on dfs_index to preserve original stable tie behavior exactly
                perimeter_table.sort(key=attrgetter('dfs_index'))
                return perimeter_table
            
            def _build_perimeter_table_acc_2(perimeter_verts: Set[TwinGraph.Vert],
                                       lining_verts: Set[TwinGraph.Vert]) -> List[BoundaryElementNT]:
                # Pre-allocate list with known size
                perimeter_table: List[BoundaryElementNT] = [None] * (len(perimeter_verts) * 2 + len(lining_verts) * 2) # type: ignore[list-item]
                be = BoundaryElementNT

                i = 0
                for v in perimeter_verts:
                    mt = v.map_tree_vert
                    perimeter_table[i] = be(mt.dfs_in, PERIMETER, v)
                    i += 1
                for v in perimeter_verts:
                    mt = v.map_tree_vert
                    perimeter_table[i] = be(mt.dfs_out, PERIMETER, v)
                    i += 1
                for v in lining_verts:
                    mt = v.map_tree_vert
                    perimeter_table[i] = be(mt.dfs_in, LINING, v)
                    i += 1
                for v in lining_verts:
                    mt = v.map_tree_vert
                    perimeter_table[i] = be(mt.dfs_out, LINING, v)
                    i += 1
                # Sort boundary elements by dfs index to group path pairs
                # Use a key on dfs_index to preserve original stable tie behavior exactly
                _perimeter_table_sort(perimeter_table)
                return perimeter_table

            def _build_perimeter_table_acc(perimeter_verts: Set[TwinGraph.Vert],
                                       lining_verts: Set[TwinGraph.Vert]) -> List[BoundaryElementNT]:
                # Build in-place with pre-sized list to avoid temporary lists and appends
                be = BoundaryElementNT
                pv_len = len(perimeter_verts)
                lv_len = len(lining_verts)
                total = (pv_len + lv_len) * 2
                perimeter_table: List[BoundaryElementNT] = [None] * total  # type: ignore[list-item]

                # Fast path: use array lookups if available, else fallback to map_tree_vert
                graph = self.region_tree.graph
                # Prefer new unified arrays; keep legacy names as fallback
                dual_in = getattr(graph, 'idx_to_dfs_in', None)
                dual_out = getattr(graph, 'idx_to_dfs_out', None)

                i = 0
                # print("Building perimeter table with", pv_len, "perimeter verts and", lv_len, "lining verts.")

                if dual_in is not None and dual_out is not None:
                    # Use index-based array lookups
                    # Preserve original build order to maintain stable tie behavior on sort(key=dfs_index):
                    # 1) all perimeter dfs_in, 2) all perimeter dfs_out, 3) all lining dfs_in, 4) all lining dfs_out
                    for v in perimeter_verts:
                        vi = v.index
                        perimeter_table[i] = be(dual_in[vi], PERIMETER, v)
                        i += 1
                    for v in perimeter_verts:
                        vi = v.index
                        perimeter_table[i] = be(dual_out[vi], PERIMETER, v)
                        i += 1
                    for v in lining_verts:
                        vi = v.index
                        perimeter_table[i] = be(dual_in[vi], LINING, v)
                        i += 1
                    for v in lining_verts:
                        vi = v.index
                        perimeter_table[i] = be(dual_out[vi], LINING, v)
                        i += 1
                else:
                    # Fallback to attribute access on map_tree_vert
                    for v in perimeter_verts:
                        mt = v.map_tree_vert
                        perimeter_table[i] = be(mt.dfs_in, PERIMETER, v)
                        i += 1
                    for v in perimeter_verts:
                        mt = v.map_tree_vert
                        perimeter_table[i] = be(mt.dfs_out, PERIMETER, v)
                        i += 1
                    for v in lining_verts:
                        mt = v.map_tree_vert
                        perimeter_table[i] = be(mt.dfs_in, LINING, v)
                        i += 1
                    for v in lining_verts:
                        mt = v.map_tree_vert
                        perimeter_table[i] = be(mt.dfs_out, LINING, v)
                        i += 1

                # Sort boundary elements by dfs index to group path pairs (dfs_index is first field)
                _perimeter_table_sort(perimeter_table)
                return perimeter_table

            def _perimeter_table_sort(perimeter_table: List[BoundaryElementNT]):
                # Use key on dfs_index to preserve exact stable tie behavior across build methods
                perimeter_table.sort(key=lambda be: be.dfs_index)

            def _collect_dual_tree_path_pairs(perimeter_table: List[BoundaryElementNT]) -> List[Tuple[BoundaryElementNT, BoundaryElementNT]]:
                # Each pair that does not have a perimeter vert inbetween is an interior path
                # We only check index-increasing pairs as we would otherwise double count
                # This means we do not need to wrap around the list
                dual_tree_path_pairs_temp: List[Tuple[BoundaryElementNT, BoundaryElementNT]] = []
                # Collect valid pairs
                for leading_idx, trailing_idx in pairwise(perimeter_table):
                    if leading_idx.boundary_type == LINING and trailing_idx.boundary_type == LINING:
                        dual_tree_path_pairs_temp.append((leading_idx, trailing_idx))
                    elif leading_idx.boundary_type == LINING and trailing_idx.boundary_type == PERIMETER:
                        dual_tree_path_pairs_temp.append((leading_idx, leading_idx)) # Pair for single lining vert abutting perimeter
                if perimeter_table[-2].boundary_type == PERIMETER and perimeter_table[-1].boundary_type == LINING:
                    dual_tree_path_pairs_temp.append((perimeter_table[-1], perimeter_table[-1])) # Pair for single lining vert abutting perimeter
                return dual_tree_path_pairs_temp

            def _consolidate_ranges_and_count_weight(dual_tree_path_pairs_temp: List[Tuple[BoundaryElementNT, BoundaryElementNT]]) -> Tuple[List[Tuple[BoundaryElementNT, BoundaryElementNT]], int]:
                dual_tree_path_pairs: List[Tuple[BoundaryElementNT, BoundaryElementNT]] = []
                double_weight = 0
                if len(dual_tree_path_pairs_temp) == 0:
                    return (dual_tree_path_pairs, double_weight) # No interior verts available

                # Consolidate adjacent ranges and count weight
                cursor_pair = dual_tree_path_pairs_temp[0]
                for pair in dual_tree_path_pairs_temp[1:]:
                    if pair[0].dfs_index == cursor_pair[1].dfs_index:
                        # Adjacent ranges, consolidate
                        cursor_pair = (cursor_pair[0], pair[1])
                    else:
                        # Not adjacent, finalize current pair
                        dual_tree_path_pairs.append(cursor_pair)
                        leading_idx, trailing_idx = cursor_pair
                        double_weight += trailing_idx.dfs_index - leading_idx.dfs_index + 1 # Add 1 to include trailing index
                        # Advance cursor
                        cursor_pair = pair
                # Finalize last pair
                dual_tree_path_pairs.append(cursor_pair)
                leading_idx, trailing_idx = cursor_pair
                double_weight += trailing_idx.dfs_index - leading_idx.dfs_index + 1 # Add 1 to include trailing index
                return (dual_tree_path_pairs, double_weight)

            def _select_checkpoint_index(double_weight: int) -> Optional[int]:
                # Select random
                if double_weight == 0:
                    return None # No interior verts available
                return random.randint(0, double_weight - 1)

            def _locate_entry_and_target_idx(dual_tree_path_pairs: List[Tuple[BoundaryElementNT, BoundaryElementNT]],
                                             checkpoint_range_target_idx: int) -> Tuple[Optional[TwinGraph.Vert], Optional[int]]:
                dfs_target_idx: Optional[int] = None
                entry_vert: Optional[TwinGraph.Vert] = None
                logical_range_covered = 0
                for leading, trailing in dual_tree_path_pairs:
                    if logical_range_covered + (trailing.dfs_index - leading.dfs_index + 1) > checkpoint_range_target_idx:
                        # Target is within this range
                        offset_into_range = checkpoint_range_target_idx - logical_range_covered
                        dfs_target_idx = leading.dfs_index + offset_into_range
                        entry_vert = leading.vert
                        break
                    logical_range_covered += trailing.dfs_index - leading.dfs_index + 1 # Add 1 to include trailing index
                return (entry_vert, dfs_target_idx)

            def _traverse_to_interior(entry_vert: TwinGraph.Vert, dfs_target_idx: int) -> TwinGraph.Vert:
                # Traverse to interior
                current_vert = entry_vert.map_tree_vert
                # print("DFS Target Index:", dfs_target_idx)
                while True:
                    # print("At dfs", current_vert.dfs_in, current_vert.dfs_out, "on vert", current_vert.twin_vert.id_str)
                    if current_vert.dfs_in == dfs_target_idx or \
                       current_vert.dfs_out == dfs_target_idx:
                        break
                    # next_vert_found = False
                    contains = current_vert.dfs_in <= dfs_target_idx <= current_vert.dfs_out
                    if contains:
                        # Target is in subtree, go deeper
                        if len(current_vert.out_edges) == 0:
                            raise ValueError("RegionTree.Region.get_uniform_interior_vert could not go down from leaf.")
                        for edge in current_vert.out_edges:
                            next_vert = edge.get_dest_from(current_vert)
                            if next_vert.dfs_in <= dfs_target_idx <= next_vert.dfs_out:
                                current_vert = next_vert
                                # print("Moving to dfs", current_vert.dfs_in, "on vert", current_vert.twin_vert.id_str)
                                break
                    else:
                        # Target is outside subtree, go up
                        # print("Targeting dfs", dfs_target_idx, "from dfs", current_vert.dfs_in, current_vert.dfs_out)
                        if current_vert.in_edge is None:
                            raise ValueError("RegionTree.Region.get_uniform_interior_vert could not go up from root.")
                        next_vert = current_vert.in_edge.get_dest_from(current_vert)
                        current_vert = next_vert
                        # print("Moving to dfs", current_vert.dfs_in, "on vert", current_vert.twin_vert.id_str)
                return current_vert.twin_vert

            # This function utilizes the dfs indexing of the dual map tree to
            #  uniformly select an interior vert.  As trees only contain one
            #  path between any two verts, if we take dfs_in and dfs_out on
            #  each lining vert, the pairs that are adjacent in sorted order 
            #  without perimeter verts inbetween define ranges within the dfs 
            #  tree that do not cross the perimeter.  As each index is a point 
            #  of crossing the perimeter, the sorted order of the indices
            #  alternates between interior and exterior paths.  This function
            #  identifies whether even or odd ranges are interior, then selects
            #  a random index within the interior ranges.
            
            lining_verts, perimeter_verts = self.get_interior_lining_verts()
            # Use accelerated builder (pre-sized + array lookups) to reduce overhead; preserves tie behavior
            perimeter_table = _build_perimeter_table(perimeter_verts, lining_verts)
            dual_tree_path_pairs_temp = _collect_dual_tree_path_pairs(perimeter_table)
            if len(dual_tree_path_pairs_temp) == 0:
                return None # No interior verts available

            dual_tree_path_pairs, double_weight = _consolidate_ranges_and_count_weight(dual_tree_path_pairs_temp)
            checkpoint_range_target_idx = _select_checkpoint_index(double_weight)
            if checkpoint_range_target_idx is None:
                return None # No interior verts available

            entry_vert, dfs_target_idx = _locate_entry_and_target_idx(dual_tree_path_pairs, checkpoint_range_target_idx)
            if dfs_target_idx is None:
                raise ValueError("RegionTree.Region.get_uniform_interior_vert could not locate dfs target index.")
            if entry_vert is None:
                raise ValueError("RegionTree.Region.get_uniform_interior_vert could not locate entry vert.")

            return _traverse_to_interior(entry_vert, dfs_target_idx)
        
        def get_uniform_interior_vert_np(self) -> Optional[TwinGraph.Vert]:
            """Hybrid NumPy implementation of uniform interior vert selection.

            Semantics match get_uniform_interior_vert exactly:
            - Builds boundary event table in four ordered blocks (perim in, perim out, lining in, lining out)
            - Stable sorting on dfs_index preserves tie ordering from build blocks
            - Interior path pair construction rules identical (lining-lining, lining-perimeter, final perimeter-lining case)
            - Consolidation of adjacent ranges where trailing_start == leading_end
            - Weight counting uses (end - start + 1) inclusive logic
            - Random checkpoint chosen uniformly in total interior index measure
            - Traversal identical: navigate map tree by dfs_in/out containment until hitting target index

            Hybrid approach: boundary table dfs values, types, and vert indices stored in NumPy arrays for fast sorting,
            but range consolidation & traversal kept in Python for clarity and minimal overhead at typical region sizes.
            Falls back to original method if NumPy dfs arrays are unavailable.
            """
            # Fallback on missing perimeter (whole graph) or absent NumPy arrays
            if len(self.dual_perimeter) == 0:
                return random.choice(list(self.region_tree.graph.dualVerts))

            graph = self.region_tree.graph
            dfs_in_np = getattr(graph, 'idx_to_dfs_in_np', None)
            dfs_out_np = getattr(graph, 'idx_to_dfs_out_np', None)
            dual_list = getattr(graph, 'dualVerts', None)
            if dfs_in_np is None or dfs_out_np is None or dual_list is None:
                return self.get_uniform_interior_vert()  # NumPy data unavailable, use legacy path

            # Tags: 0 = PERIMETER, 1 = LINING (match original int tags)
            PERIMETER = 0
            LINING = 1

            lining_verts, perimeter_verts = self.get_interior_lining_verts()
            if len(lining_verts) == 0:
                return None  # No interior verts available

            # Build index arrays for perimeter and lining vert indices
            # Note: Order of iteration over sets matches original (arbitrary set order) -> preserved via stable sort
            pv_idx = np.fromiter((v.index for v in perimeter_verts), dtype=np.int32, count=len(perimeter_verts))
            lv_idx = np.fromiter((v.index for v in lining_verts), dtype=np.int32, count=len(lining_verts))

            # Boundary dfs indices in block order matching original builder
            boundary_dfs = np.concatenate([
                dfs_in_np[pv_idx],
                dfs_out_np[pv_idx],
                dfs_in_np[lv_idx],
                dfs_out_np[lv_idx]
            ])

            # Boundary types aligned with blocks (stable ordering preserved by mergesort argsort)
            boundary_type = np.concatenate([
                np.zeros(len(pv_idx), dtype=np.int8),  # perim in
                np.zeros(len(pv_idx), dtype=np.int8),  # perim out
                np.ones(len(lv_idx), dtype=np.int8),   # lining in
                np.ones(len(lv_idx), dtype=np.int8)    # lining out
            ])

            # Vert indices (two events per vertex: in/out)
            vert_indices = np.concatenate([
                pv_idx,
                pv_idx,
                lv_idx,
                lv_idx
            ])

            # Stable sort by dfs index to group path crossings; mergesort guarantees stability
            order = np.argsort(boundary_dfs, kind='mergesort')
            boundary_dfs = boundary_dfs[order]
            boundary_type = boundary_type[order]
            vert_indices = vert_indices[order]

            n = boundary_dfs.size
            if n == 0:
                return None

            # Vectorized detection of candidate pairs (interior path intervals)
            # Pair types:
            #   lining-lining -> (start_i, end_{i+1})
            #   lining-perimeter -> single vertex interval (start_i, start_i)
            # Final special case: perimeter-lining at end -> single vertex interval
            leading_type = boundary_type[:-1]
            trailing_type = boundary_type[1:]
            leading_dfs = boundary_dfs[:-1]
            trailing_dfs = boundary_dfs[1:]

            lining_lining_mask = (leading_type == LINING) & (trailing_type == LINING)
            lining_perim_mask = (leading_type == LINING) & (trailing_type == PERIMETER)

            # Collect starts/ends and positions of leading elements
            pair_starts = []  # dfs index of leading event
            pair_ends = []    # dfs index of trailing event (or same for single)
            leading_positions = []  # array index of leading event (for entry vert retrieval)

            # Lining-lining intervals
            ll_indices = np.nonzero(lining_lining_mask)[0]
            for i in ll_indices:
                pair_starts.append(leading_dfs[i])
                pair_ends.append(trailing_dfs[i])
                leading_positions.append(i)

            # Lining-perimeter singletons
            lp_indices = np.nonzero(lining_perim_mask)[0]
            for i in lp_indices:
                pair_starts.append(leading_dfs[i])
                pair_ends.append(leading_dfs[i])  # single vertex interval
                leading_positions.append(i)

            # Final perimeter-lining singleton case
            if n >= 2 and boundary_type[-2] == PERIMETER and boundary_type[-1] == LINING:
                pair_starts.append(boundary_dfs[-1])
                pair_ends.append(boundary_dfs[-1])
                leading_positions.append(n - 1)

            if len(pair_starts) == 0:
                return None  # No interior verts available

            # Consolidate adjacent intervals where start of next equals end of current
            consolidated = []  # list of (start_dfs, end_dfs, leading_pos)
            cursor_start = pair_starts[0]
            cursor_end = pair_ends[0]
            cursor_lead_pos = leading_positions[0]
            double_weight = 0
            for s, e, lp in zip(pair_starts[1:], pair_ends[1:], leading_positions[1:]):
                if s == cursor_end:  # adjacency -> extend range
                    cursor_end = e
                else:
                    consolidated.append((cursor_start, cursor_end, cursor_lead_pos))
                    double_weight += (cursor_end - cursor_start + 1)
                    cursor_start, cursor_end, cursor_lead_pos = s, e, lp
            # Final interval
            consolidated.append((cursor_start, cursor_end, cursor_lead_pos))
            double_weight += (cursor_end - cursor_start + 1)

            if double_weight == 0:
                return None

            # Select checkpoint uniformly across inclusive dfs index measure
            checkpoint = random.randint(0, double_weight - 1)

            # Locate interval containing checkpoint
            logical_covered = 0
            entry_vert_index = None
            dfs_target = None
            for start_dfs, end_dfs, lead_pos in consolidated:
                span = end_dfs - start_dfs + 1
                if logical_covered + span > checkpoint:
                    offset = checkpoint - logical_covered
                    dfs_target = start_dfs + offset
                    entry_vert_index = vert_indices[lead_pos]
                    break
                logical_covered += span

            if dfs_target is None or entry_vert_index is None:
                raise ValueError("RegionTree.Region.get_uniform_interior_vert_np could not locate dfs target or entry vert.")

            # dual_list is a set (unordered, non-subscriptable). Provide cached index->vertex map.
            # TwinGraph uses __slots__, so we cannot attach a new attribute to cache.
            # Build a local mapping (cost negligible compared to perimeter table construction).
            vert_lookup = {}
            for v in dual_list:  # unique indices assumed
                vert_lookup[v.index] = v
            entry_vert = vert_lookup.get(entry_vert_index)
            if entry_vert is None:
                raise ValueError(
                    f"RegionTree.Region.get_uniform_interior_vert_np could not resolve vertex index {entry_vert_index} among {len(vert_lookup)} dual verts."
                )

            # Traversal identical to original logic
            def _traverse(entry_vert: TwinGraph.Vert, dfs_target_idx: int) -> TwinGraph.Vert:
                current_vert = entry_vert.map_tree_vert
                while True:
                    if current_vert.dfs_in == dfs_target_idx or current_vert.dfs_out == dfs_target_idx:
                        break
                    contains = current_vert.dfs_in <= dfs_target_idx <= current_vert.dfs_out
                    if contains:
                        if len(current_vert.out_edges) == 0:
                            raise ValueError("RegionTree.Region.get_uniform_interior_vert_np could not go down from leaf.")
                        for edge in current_vert.out_edges:
                            next_vert = edge.get_dest_from(current_vert)
                            if next_vert.dfs_in <= dfs_target_idx <= next_vert.dfs_out:
                                current_vert = next_vert
                                break
                    else:
                        if current_vert.in_edge is None:
                            raise ValueError("RegionTree.Region.get_uniform_interior_vert_np could not go up from root.")
                        current_vert = current_vert.in_edge.get_dest_from(current_vert)
                return current_vert.twin_vert

            return _traverse(entry_vert, dfs_target)

        def get_central_region_vert(self) -> Optional[TwinGraph.Vert]:
            """
            Conducts a breadth-first search from the dual perimeter to find a central vert.
            Returns None if no interior vert is found (region is a single primal vert).
            Assumes that the region is connected.
            """
            visited = self.get_perimeter_verts()
            lining_verts, _ = self.get_interior_lining_verts()
            queue = deque(lining_verts)
            dual_vert = None
            while queue:
                dual_vert = queue.popleft()
                if dual_vert in visited:
                    continue
                visited.add(dual_vert)
                for neighbor_edge in dual_vert.cc_edges:
                    neighbor, _ = neighbor_edge.get_dual_dest_from(dual_vert)
                    if neighbor not in visited:
                        queue.append(neighbor)
            return dual_vert

        def get_walk_random_interior_vert(self, steps: int = 100) -> Optional[TwinGraph.Vert]:
            # Special case: no perimeter means whole-graph; perform an unrestricted walk
            if len(self.dual_perimeter) == 0:
                vert = random.choice(list(self.region_tree.graph.dualVerts))
                return vert

            # Build perimeter vertex set to avoid stepping onto boundary
            perimeter_verts = self.get_perimeter_verts()

            # Helpers ------------------------------------------------------
            def _find_interior_start_vert() -> Optional[TwinGraph.Vert]:
                """Find an interior start vertex adjacent to the perimeter without landing on a perimeter vert."""
                start_idx = random.randint(0, len(self.dual_perimeter) - 1)
                offset = 0
                while offset < len(self.dual_perimeter):
                    boundary_start_edge, boundary_start_dir = self.dual_perimeter[(start_idx + offset) % len(self.dual_perimeter)]
                    boundary_next_edge, boundary_next_dir = self.dual_perimeter[(start_idx + offset + 1) % len(self.dual_perimeter)]
                    _, boundary_next_vert = boundary_next_edge.get_dual_vert_pair(boundary_next_dir)
                    _, rotary_center = boundary_start_edge.get_dual_vert_pair(boundary_start_dir)
                    
                    if rotary_center is None:
                        raise ValueError("RegionTree.Region.get_walk_random_interior_vert found edge with undefined dual vert.")

                    while True:
                        next_edge, _ = boundary_start_edge.get_dual_cc_next_edge(rotary_center)
                        candidate_interior_vert, _ = next_edge.get_dual_dest_from(rotary_center)

                        # Vert in sweep that lies between boundary_start and boundary_next is interior
                        if candidate_interior_vert is not boundary_next_vert:
                            return candidate_interior_vert

                        if next_edge == boundary_next_edge:
                            break

                    offset += 1
                return None

            def _walk_interior_optimistic(start: TwinGraph.Vert, step_cap: int) -> TwinGraph.Vert:
                """
                Optimistic Walk: Blindly picks a direction. 
                If it hits a wall (perimeter), it 'rewinds' (stays put) and burns the step.
                """
                current = start
                for _ in range(max(0, step_cap)):
                    edges = current.cc_edges
                    degree = len(edges)
                    if degree == 0:
                        break

                    # 1. Blindly pick a random edge
                    # random.choice is slightly faster than randint + indexing for lists in Python
                    edge = random.choice(edges) 
                    
                    # 2. Step
                    candidate, _ = edge.get_dual_dest_from(current)

                    # 3. Check (Rewind if invalid)
                    if candidate not in perimeter_verts:
                        current = candidate
                    # else: 
                    #   Hit a wall. Effectively "rewind" by keeping 'current' as is.
                    #   We accept the 'wasted' step to avoid the cost of finding a valid one.
                
                return current

            # Execution ------------------------------------------------------
            
            src_vert: Optional[TwinGraph.Vert] = _find_interior_start_vert()

            if src_vert is None:
                return None

            current = _walk_interior_optimistic(src_vert, steps)

            return current
        
        def check_center(self) -> bool:
            is_center = True

            for edge in self.region_edges:
                if edge.ab_weight_differential is None:
                    raise ValueError("RegionTree.Region.check_center found edge with undefined ab_weight_differential.")
                
                opposite_weight = edge.get_opposite_side_weight(self)
                midpoint_weight = self.region_tree.total_weight / 2
                epsilon_buffer = self.region_tree.epsilon * midpoint_weight
                if opposite_weight is None:
                    raise ValueError("RegionTree.Region.check_center found edge with undefined opposite side weight.")
                elif opposite_weight >= midpoint_weight + epsilon_buffer:
                    is_center = False
                # print(f"Edge {edge.id_str} opposite weight {opposite_weight}, midpoint {midpoint_weight}, region weight {self.weight}")

            if is_center:
                if self.validate_split_possible():
                    # TODO: Remove profiling
                    # old_region_size = self.region_tree.central_region.weight if self.region_tree.central_region is not None else None
                    # if old_region_size is not None:
                    #     new_region_size = self.weight
                        # print(f"Ratio: {old_region_size} -> {new_region_size} = {new_region_size / old_region_size:.4f}")

                    # self.region_tree.target_regions.append(self)
                    return True
                    # print("Region", self.id_str, "is central region of region tree", self.region_tree.id_str, "with", self.weight, "weight.")
                # else:
                #     self.region_tree.central_regions = None
                    # print("Region", self.id_str, "is central, but a 2-split not possible.")
            return False

        def validate_split_possible(self) -> bool:
            if len(self.cc_region_edges) == 0:
                return True # No edges implies whole graph, which is splittable

            src_graph_weight = self.region_tree.total_weight
            midpoint_weight = src_graph_weight / 2
            epsilon_buffer = self.region_tree.epsilon * midpoint_weight

            trailing_idx = 0 # Inclusive
            leading_idx = 1 # Exclusive
            a_1_weight = self.cc_region_edges[0].get_opposite_side_weight(self)
            if a_1_weight is None:
                raise ValueError("RegionTree.Region.validate_split_possible found edge with undefined opposite side weight.")
            a_2_weight = src_graph_weight - self.weight - a_1_weight
            
            while True:
                if a_1_weight <= midpoint_weight + epsilon_buffer and a_2_weight <= midpoint_weight + epsilon_buffer:
                    return True
                # Terminate when trailing has wrapped around or when leading has caught trailing
                if trailing_idx == len(self.cc_region_edges) or (leading_idx > len(self.cc_region_edges) and leading_idx % len(self.cc_region_edges) == trailing_idx):
                    return False
                
                if a_1_weight > midpoint_weight + epsilon_buffer:
                    # Move trailing forward
                    dropped_region_weight = self.cc_region_edges[trailing_idx % len(self.cc_region_edges)].get_opposite_side_weight(self)
                    if dropped_region_weight is None:
                        raise ValueError("RegionTree.Region.validate_split_possible found edge with undefined opposite side weight.")
                    a_1_weight -= dropped_region_weight
                    a_2_weight += dropped_region_weight
                    trailing_idx += 1
                else:
                    # Move leading forward
                    added_region_weight = self.cc_region_edges[leading_idx % len(self.cc_region_edges)].get_opposite_side_weight(self)
                    if added_region_weight is None:
                        raise ValueError("RegionTree.Region.validate_split_possible found edge with undefined opposite side weight.")
                    a_1_weight += added_region_weight
                    a_2_weight -= added_region_weight
                    leading_idx += 1

    class Edge:
        __slots__ = (
            'end_A', 
            'end_B', 
            'twin_graph_edge', 
            'ab_weight_differential', 
            'region_tree', 
            'id_str'
        )

        # Edge A and B ends must match directionally with underlying twin_graph_edge
        end_A: RegionTree.Region
        end_B: RegionTree.Region
        twin_graph_edge: TwinGraph.QuadEdge
        ab_weight_differential: Optional[int] # Difference in weight between all regions on A side and all regions on B side (weight of A - weight of B)

        # Tree
        region_tree: RegionTree # Must be set when edge is added to tree

        # Instance labeling
        instance_counter: int = 0
        id_str: str

        def __init__(self) -> None:
            self.id_str = f'RE{RegionTree.Edge.instance_counter:x}' # Hexadecimal ID for easier reading in debug
            RegionTree.Edge.instance_counter += 1

        # Get the vertex on the other end of an edge along with directional annotation
        def get_dest_from(self, src: RegionTree.Region) -> Tuple[RegionTree.Region, RegionTree.EdgeDir]:
            if src == self.end_A:
                return (self.end_B, TwinGraph.EdgeDir.AB)
            elif src == self.end_B:
                return (self.end_A, TwinGraph.EdgeDir.BA)
            else:
                raise KeyError("Src for get_primal_dest_from is not on edge.")

        # Get the angle to the vertex on the other end of an edge along with directional annotation
        # def get_rad_from(self, src: RegionTree.Region) -> Tuple[float, RegionTree.EdgeDir]:
        #     _, dir = self.get_dest_from(src)
        #     src_rad = self.twin_graph_edge.get_dual_rad_along(dir)

        #     return (src_rad, dir)

        # Get the total weight on the opposite side of the edge of the given region
        # Returns None if ab_weight_differential is not yet calculated
        def get_opposite_side_weight(self, src: RegionTree.Region) -> Optional[int]:
            if self.ab_weight_differential is None:
                return None
            
            assert (self.ab_weight_differential + self.region_tree.graph.total_primal_weight) % 2 == 0, "RegionTree.Edge ab_weight_differential and graph total_primal_weight do not sum to an even number, cannot divide by 2 cleanly to get total weight a."
            assert (self.region_tree.graph.total_primal_weight - self.ab_weight_differential) % 2 == 0, "RegionTree.Edge region_tree.graph total_primal_weight and ab_weight_differential do not subtract to an even number, cannot divide by 2 cleanly to get total weight b."

            total_weight_a = (self.ab_weight_differential + self.region_tree.graph.total_primal_weight) // 2
            total_weight_b = (self.region_tree.graph.total_primal_weight - self.ab_weight_differential) // 2

            if src == self.end_A:
                return total_weight_b
            elif src == self.end_B:
                return total_weight_a
            else:
                raise KeyError("Src for get_opposite_side_weight is not on edge.")

        # Calculate the weight differential based on the side of the graph pointed to by dir
        def calculate_weight_differential_from_dir(self, dir: TwinGraph.EdgeDir) -> None:
            if dir == TwinGraph.EdgeDir.BA:
                side_A_weight = self.end_A.weight
                for edge in self.end_A.region_edges - {self}:
                    weight = edge.get_opposite_side_weight(self.end_A)
                    if weight is not None:
                        side_A_weight += weight
                    else:
                        raise ValueError("Cannot calculate weight differential, side A of edge has undefined weight.")
                side_B_weight = self.region_tree.graph.total_primal_weight - side_A_weight
            elif dir == TwinGraph.EdgeDir.AB:
                side_B_weight = self.end_B.weight
                for edge in self.end_B.region_edges - {self}:
                    weight = edge.get_opposite_side_weight(self.end_B)
                    if weight is not None:
                        side_B_weight += weight
                    else:
                        raise ValueError("Cannot calculate weight differential, side B of edge has undefined weight.")
                side_A_weight = self.region_tree.graph.total_primal_weight - side_B_weight

            self.ab_weight_differential = side_A_weight - side_B_weight

            # Check for edge center
            if abs(self.ab_weight_differential) <= self.region_tree.epsilon * self.region_tree.total_weight: # 2 * self.region_tree.epsilon * (1/2) * self.region_tree.total_weight
                # if len(self.region_tree.edge_centers) > 0:
                #     raise ValueError("Cannot set edge as edge center, region tree already has edge center.")

                self.region_tree.edge_centers.add(self)
                # self.region_tree.central_region = None
                # print("Edge", self.id_str, "is edge center of region tree", self.region_tree.id_str)


            
