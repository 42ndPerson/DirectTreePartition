from __future__ import annotations
import bisect
import random
from typing import List, Tuple, Set
from enum import Enum
import warnings

from TwinGraph import *
from Euclid import *

class RegionTree:
    __slots__ = (
        'regions', 
        'edges', 
        'graph', 
        'central_region', 
        'edge_center',  
        'id_str'
    )

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
        total_weight = graph.total_primal_weight
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
            verts = set()
            for edge, dir in self.dual_perimeter:
                src, _ = edge.get_dual_vert_pair(dir)
                verts.add(src)
            return verts
        
        # def get_perimeter_edges(self) -> Set[TwinGraph.QuadEdge]:
        #     return {edge for edge, _ in self.dual_perimeter}

        def get_interior_lining_verts(self) -> Set[TwinGraph.Vert]:
            if len(self.dual_perimeter) == 0:
                return {e.get_dual_dest_from(self.region_tree.graph.external_dual_vert)[0] 
                        for e in self.region_tree.graph.external_dual_vert.cc_edges}

            full_perimeter_edges = self.dual_perimeter_edges #self.get_perimeter_edges()
            lining_verts: Set[TwinGraph.Vert] = set()

            for edge, dir in self.dual_perimeter:
                rotary_center: TwinGraph.Vert
                _, rotary_center = edge.get_dual_vert_pair(dir)

                working_edge = edge
                while True:
                    working_edge, _ = working_edge.get_dual_cc_next_edge(rotary_center)
                    if working_edge in full_perimeter_edges: 
                        break

                    lining_verts.add(working_edge.get_dual_dest_from(rotary_center)[0])
                    # Performance note: Set is absorbing ~50% double vert detections with grid graphs
                    
            return lining_verts

        def get_uniform_interior_vert(self) -> Optional[TwinGraph.Vert]:
            if len(self.dual_perimeter) == 0:
                return random.choice(list(self.region_tree.graph.dualVerts))

            full_perimeter_verts = [e.get_dual_vert_pair(dir)[0] for e, dir in self.dual_perimeter]
            print("Perimeter Edges:", [e.id_str for e, _ in self.dual_perimeter])
            print("Perim Verts:", [v.id_str for v in full_perimeter_verts])
            print("Measured weight:", self.region_tree.graph.count_primal_verts_within_perim(self.dual_perimeter))
            dual_tree_checkpoints: List[int] = []
            dual_tree_checkpoints_interior_flowing: List[bool] = []
            dual_tree_checkpoints_map: Dict[int, TwinGraph.Vert] = {}
            dual_tree_path_pairs: List[Tuple[int, int]] = []

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

            # Collect all entry and exit points to region
            for edge, dir in self.dual_perimeter:
                _, vert = edge.get_dual_vert_pair(dir)

                dual_tree_checkpoints.append(vert.map_tree_vert.dfs_in)
                dual_tree_checkpoints.append(vert.map_tree_vert.dfs_out)
                dual_tree_checkpoints_map[vert.map_tree_vert.dfs_in] = vert
                dual_tree_checkpoints_map[vert.map_tree_vert.dfs_out] = vert

                # Determine whether interior or exterior
                # Start by extracting interior-side dual-tree edges
                interior_dual_mapping_tree_edges: Set[TwinGraph.QuadEdge] = set()
                rotary_center = vert
                current_edge = edge
                while True:
                    # Advance sweep
                    current_edge, _ = current_edge.get_dual_cc_next_edge(rotary_center)

                    # Check whether sweep has re-encountered perimeter edge
                    if current_edge in self.dual_perimeter_edges:
                        break

                    # Check whether edge is in mapping tree
                    if current_edge in self.region_tree.graph.dualMap.quad_edges:
                        interior_dual_mapping_tree_edges.add(current_edge)
                # Determine flow

            
            print("Raw checkpoints:", dual_tree_checkpoints)
            print("Raw interior flowing:", dual_tree_checkpoints_interior_flowing)
            # Sort to identify index ranges for paths between checkpoints
            zipped_checkpoint_info = list(zip(dual_tree_checkpoints, dual_tree_checkpoints_interior_flowing))
            print("Zipped checkpoints:", zipped_checkpoint_info)
            zipped_checkpoint_info.sort() # Group subranges that identify paths between perimeter checkpoints
            dual_tree_checkpoints, dual_tree_checkpoints_interior_flowing = map(list, zip(*zipped_checkpoint_info))
            print("Sorted checkpoints:", dual_tree_checkpoints)
            print("Sorted interior flowing:", dual_tree_checkpoints_interior_flowing)

            # Count interior dual verts
            print(dual_tree_checkpoints)
            double_weight = 0
            for perim_idx in range(len(dual_tree_checkpoints)):
                next_perim_idx = (perim_idx + 1) % len(dual_tree_checkpoints)
                if dual_tree_checkpoints_interior_flowing[perim_idx]:
                    double_weight += dual_tree_checkpoints[next_perim_idx] - dual_tree_checkpoints[perim_idx] - 1 # Subtract 1 to not count the perimeter verts

            # Select random
            print("Double Weight:", double_weight)
            if double_weight == 0:
                return None # No interior verts available
            checkpoint_range_target_idx = random.randint(0, double_weight - 1)
            dfs_target_idx = None

            # Locate perimeter entry point
            entry_vert: Optional[TwinGraph.Vert] = None
            logical_range_covered = 0
            print(dual_tree_checkpoints)
            for perim_idx in range(0, len(dual_tree_checkpoints), 2):
                if dual_tree_checkpoints_interior_flowing[perim_idx]:
                    range_start_checkpoint_dfs_idx = dual_tree_checkpoints[perim_idx]
                    range_end_checkpoint_dfs_idx = dual_tree_checkpoints[(perim_idx + 1) % len(dual_tree_checkpoints)]
                    range_size = range_end_checkpoint_dfs_idx - range_start_checkpoint_dfs_idx - 1 # Subtract 1 to not count the perimeter verts

                    # Check if target index is within this range
                    if checkpoint_range_target_idx >= logical_range_covered and \
                    checkpoint_range_target_idx < logical_range_covered + range_size:
                        entry_vert = dual_tree_checkpoints_map[range_start_checkpoint_dfs_idx]

                        # Get dfs index of target interior vert
                        dfs_target_idx = \
                            range_start_checkpoint_dfs_idx + 1 + (checkpoint_range_target_idx - logical_range_covered)
                            # Add 1 to skip perimeter vert

                        break
                    
                    # Track logical range covered
                    logical_range_covered += range_size

            # Traverse to interior
            print("Target dfs idx:", dfs_target_idx)
            print("Entry vert:", entry_vert.id_str if entry_vert is not None else "None")
            if entry_vert is None:
                raise ValueError("RegionTree.Region.get_uniform_interior_vert could not locate entry vert.")
            current_vert = entry_vert
            while True:
                if current_vert.map_tree_vert.dfs_in == dfs_target_idx or \
                   current_vert.map_tree_vert.dfs_out == dfs_target_idx:
                    break
                
                next_vert_found = False
                for edge in current_vert.cc_edges: # Sweep through edges to find direction including target
                    next_vert, _ = edge.get_dual_dest_from(current_vert)
                    if next_vert.map_tree_vert.dfs_in <= dfs_target_idx and \
                        next_vert.map_tree_vert.dfs_out >= dfs_target_idx:
                            current_vert = next_vert
                            next_vert_found = True
                            break
                print("Stepped to vert:", current_vert.id_str)
                if not next_vert_found:
                    raise ValueError("RegionTree.Region.get_uniform_interior_vert could not find interior target along implied path.")

            return current_vert

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
                    # print("Region", self.id_str, "is central region of region tree", self.region_tree.id_str, "with", self.weight, "weight.")
                else:
                    self.region_tree.central_region = None
                    # print("Region", self.id_str, "is central, but a 2-split not possible.")

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
            if self.ab_weight_differential == 0:
                if self.region_tree.edge_center is not None:
                    raise ValueError("Cannot set edge as edge center, region tree already has edge center.")

                self.region_tree.edge_center = self
                # self.region_tree.central_region = None
                # print("Edge", self.id_str, "is edge center of region tree", self.region_tree.id_str)


            
