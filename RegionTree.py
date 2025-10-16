from __future__ import annotations
from typing import List, Tuple, Set
from enum import Enum
import warnings

from TwinGraph import *
from Euclid import *

class RegionTree:
    regions: Set[RegionTree.Region]
    edges: Set[RegionTree.Edge]

    graph: TwinGraph

    # Central region and edge
    central_region: Optional[RegionTree.Region]
    edge_center: Optional[RegionTree.Edge]

    # Instance labeling
    instance_counter: int = 0
    id_str: str

    type EdgeDir = TwinGraph.EdgeDir

    def __init__(self, graph: TwinGraph) -> None:
        # The src region of a region tree is perimter-less
        # The only perimeter it could have would be the all
        #  the edges to the dual exterior listed twice. 
        # However, this would immediately cut off any potential
        #  primal trees that connect along the edges of the 
        #  graph.
        # We are able to omit the perimeter to solve this
        #  problem because the graph is finite, so all paths
        #  are contained naturally by the graph structure.
        total_weight = sum(vert.weight for vert in graph.primalVerts)
        self.regions = set()
        self.add_region(RegionTree.Region(total_weight, [])) # Start with single region containing whole graph
        self.edges = set()

        self.graph = graph

        self.central_region = next(iter(self.regions))
        self.edge_center = None

        self.id_str = f"RTree{RegionTree.instance_counter}"
        RegionTree.instance_counter += 1

    def add_region(self, vert: RegionTree.Region) -> None:
        self.regions.add(vert) 
        vert.region_tree = self

    def remove_region(self, vert: RegionTree.Region) -> None:
        for edge in list(vert.region_edges):
            self.remove_edge(edge)
        self.regions.remove(vert)

        if vert == self.central_region:
            self.central_region = None
            print("Central region removed, no central region now set.")

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
        # Attributes
        weight: int # Total weight of the region
        dual_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] # Edges that define the perimeter of the region in the dual graph

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
            verts = set()
            for edge, dir in self.dual_perimeter:
                src, _ = edge.get_dual_vert_pair(dir)
                verts.add(src)
            return verts
        
        def get_perimeter_edges(self) -> Set[TwinGraph.QuadEdge]:
            return {edge for edge, _ in self.dual_perimeter}

        def get_interior_lining_verts(self) -> Set[TwinGraph.Vert]:
            if len(self.dual_perimeter) == 0:
                return self.region_tree.graph.dualVerts - {v for v in self.region_tree.graph.dualVerts if v.role == TwinGraph.VertRole.DUAL_EXTERIOR}

            full_perimeter_edges = self.get_perimeter_edges()
            lining_verts: Set[TwinGraph.Vert] = set()

            for edge, dir in self.dual_perimeter:
                rotary_center: TwinGraph.Vert
                if dir == TwinGraph.EdgeDir.AB and edge.dual_BA is not None:
                    rotary_center = edge.dual_BA
                elif dir == TwinGraph.EdgeDir.BA and edge.dual_AB is not None:
                    rotary_center = edge.dual_AB
                else:
                    raise ValueError("RegionTree.Region.get_interior_lining_verts found edge with invalid direction or unassigned end dual verts in dual_perimeter.")
                
                working_edge = edge
                while True:
                    working_edge, _ = working_edge.get_dual_cc_next_edge(rotary_center)
                    if working_edge in full_perimeter_edges:
                        break
                    lining_verts.add(working_edge.get_dual_dest_from(rotary_center)[0])
                        
            return lining_verts

        def check_center(self):
            is_center = True

            for edge in self.region_edges:
                if edge.ab_weight_differential is None:
                    raise ValueError("RegionTree.Region.check_center found edge with undefined ab_weight_differential.")
                
                opposite_weight = edge.get_opposite_side_weight(self)
                midpoint_weight = self.region_tree.graph.total_primal_weight / 2
                if opposite_weight is None:
                    raise ValueError("RegionTree.Region.check_center found edge with undefined opposite side weight.")
                elif opposite_weight >= midpoint_weight:
                    is_center = False
                # print(f"Edge {edge.id_str} opposite weight {opposite_weight}, midpoint {midpoint_weight}, region weight {self.weight}")

            if is_center:
                if self.validate_split_possible():
                    self.region_tree.central_region = self
                    print("Region", self.id_str, "is central region of region tree", self.region_tree.id_str, "with", self.weight, "weight.")
                else:
                    self.region_tree.central_region = None
                    print("Region", self.id_str, "is central, but a 2-split not possible.")

        def validate_split_possible(self) -> bool:
            if len(self.cc_region_edges) == 0:
                return True # No edges implies whole graph, which is splittable

            graph_weight = self.region_tree.graph.total_primal_weight
            midpoint_weight = graph_weight / 2

            trailing_idx = 0 # Inclusive
            leading_idx = 1 # Exclusive
            a_1_weight = self.cc_region_edges[0].get_opposite_side_weight(self)
            if a_1_weight is None:
                raise ValueError("RegionTree.Region.validate_split_possible found edge with undefined opposite side weight.")
            a_2_weight = graph_weight - self.weight - a_1_weight
            
            while True:
                if a_1_weight <= midpoint_weight and a_2_weight <= midpoint_weight:
                    return True
                # Terminate when trailing has wrapped around or when leading has caught trailing
                if trailing_idx == len(self.cc_region_edges) or (leading_idx > len(self.cc_region_edges) and leading_idx % len(self.cc_region_edges) == trailing_idx):
                    return False
                
                if a_1_weight > midpoint_weight:
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
        def get_rad_from(self, src: RegionTree.Region) -> Tuple[float, RegionTree.EdgeDir]:
            _, dir = self.get_dest_from(src)
            src_rad = self.twin_graph_edge.get_dual_rad_along(dir)

            return (src_rad, dir)

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
            if self.ab_weight_differential == 0:
                if self.region_tree.edge_center is not None:
                    raise ValueError("Cannot set edge as edge center, region tree already has edge center.")

                self.region_tree.edge_center = self
                print("Edge", self.id_str, "is edge center of region tree", self.region_tree.id_str)


            
