from __future__ import annotations
from collections import deque
from typing import Dict, List, Optional, Set, Tuple
import random

from TwinGraph import TwinGraph

class ReCom:
    graph: TwinGraph
    districts: Set[ReCom.District]
    adjacencies: Set[ReCom.Adjacency]

    def __init__(self, graph: TwinGraph, district_groupings: List[Set[TwinGraph.Vert]]):
        self.graph = graph
        self.districts = set()
        self.adjacencies = set()

        adjacency_lookup: Dict[TwinGraph.QuadEdge, ReCom.District] = dict()
        for grouping in district_groupings:
            adjacent_districts: Set[ReCom.District] = set()

            # Collect boundary edges
            exterior_edges: Set[TwinGraph.QuadEdge] = set()
            starting_edge: Optional[TwinGraph.QuadEdge] = None
            starting_dir: Optional[TwinGraph.EdgeDir] = None
            
            for v in grouping:
                for edge in v.cc_edges:
                    neighbor, _ = edge.get_primal_dest_from(v)
                    if neighbor not in grouping:
                        exterior_edges.add(edge)
                        if starting_edge is None:
                            starting_edge = edge
                            # Determine orientation: We want v (In) on the Right.
                            if v == edge.primal_A:
                                starting_dir = TwinGraph.EdgeDir.BA
                            else:
                                starting_dir = TwinGraph.EdgeDir.AB

            # Trace dual perimeter
            perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] = []
            
            if starting_edge is not None:
                curr_edge = starting_edge
                curr_dir = starting_dir
                
                # Start vertex of the dual edge
                assert curr_dir is not None
                curr_vert, _ = curr_edge.get_dual_vert_pair(curr_dir)
                
                # Check adjacency for start edge
                if curr_edge in adjacency_lookup:
                    adjacent_districts.add(adjacency_lookup[curr_edge])

                start_edge_fixed = curr_edge
                start_dir_fixed = curr_dir
                
                while True:
                    perimeter.append((curr_edge, curr_dir))
                    
                    # Traverse edge to next vertex
                    curr_vert, _ = curr_edge.get_dual_dest_from(curr_vert)
                    
                    # Scan CW for next boundary edge
                    scan_edge = curr_edge
                    while True:
                        scan_edge, scan_dir = scan_edge.get_dual_cc_prev_edge(curr_vert)
                        if scan_edge in exterior_edges:
                            curr_edge = scan_edge
                            curr_dir = scan_dir
                            break
                        if scan_edge == curr_edge:
                             raise Exception("Infinite loop looking for next boundary edge")

                    # Check adjacency
                    if curr_edge in adjacency_lookup:
                        adjacent_districts.add(adjacency_lookup[curr_edge])
                    
                    if curr_edge == start_edge_fixed and curr_dir == start_dir_fixed:
                        break

            # Create district
            population = self.graph.count_primal_verts_within_perim(perimeter)
            district = ReCom.District(perimeter, population)
            self.districts.add(district)
            
            # Register adjacencies
            for neighbor in adjacent_districts:
                adj = ReCom.Adjacency(district, neighbor)
                self.adjacencies.add(adj)
                district.adjacencies.add(adj)
                neighbor.adjacencies.add(adj)
            
            # Update lookup
            for edge, _ in perimeter:
                adjacency_lookup[edge] = district

    def markov_step(self) -> None:
        # Pick random adjacency
        adjacency = random.choice(list(self.adjacencies))

        # Combined region
        district_A_edges: Set[TwinGraph.QuadEdge] = set()
        district_A_perim_verts: Set[TwinGraph.Vert] = set()
        for edge, _ in adjacency.district_A.perimeter:
            district_A_edges.add(edge)
            v1, v2 = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
            district_A_perim_verts.add(v1)
            district_A_perim_verts.add(v2)
        district_B_edges: Set[TwinGraph.QuadEdge] = set()
        district_B_perim_verts: Set[TwinGraph.Vert] = set()
        for edge, _ in adjacency.district_B.perimeter:
            district_B_edges.add(edge)
            v1, v2 = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
            district_B_perim_verts.add(v1)
            district_B_perim_verts.add(v2)

        pass

    class District:
        perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]
        edge_set: Set[TwinGraph.QuadEdge]
        population: int
        adjacencies: Set[ReCom.Adjacency]

        def __init__(self, perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]], population: int):
            self.perimeter = perimeter
            self.edge_set = set(edge for edge, _ in perimeter)
            self.population = population
            self.adjacencies = set()

    class Adjacency:
        district_A: ReCom.District
        district_B: ReCom.District

        def __init__(self, district_A: ReCom.District, district_B: ReCom.District):
            self.district_A = district_A
            self.district_B = district_B