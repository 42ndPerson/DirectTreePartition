from __future__ import annotations
from typing import List, Tuple, Set

from TwinGraph import *
from Euclid import *

class RegionTree:
    regions: Set[RegionTree.Region]
    edges: Set[RegionTree.Edge]

    graph: TwinGraph

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
        self.regions = { RegionTree.Region(total_weight, []) }
        self.edges = set()

        self.graph = graph

        self.id_str = f"RTree{RegionTree.instance_counter}"
        RegionTree.instance_counter += 1

    def add_region(self, vert: RegionTree.Region) -> None: 
        print("Adding region", vert.id_str, "-----") 
        self.regions.add(vert)  

    def remove_region(self, vert: RegionTree.Region) -> None:
        print("Removing region", vert.id_str, "-----")
        for edge in list(vert.region_edges):
            self.remove_edge(edge)
        self.regions.remove(vert)

    def add_edge(self, edge: RegionTree.Edge) -> None:
        print("Adding edge", edge.id_str, "with twin", edge.twin_graph_edge.id_str)
        self.edges.add(edge)
        edge.end_A.register_edge(edge)
        edge.end_B.register_edge(edge)

    def remove_edge(self, edge: RegionTree.Edge) -> None:
        print("Removing edge", edge.id_str, "with twin", edge.twin_graph_edge.id_str)
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
        weight: float # Total weight of the region
        dual_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]] # Edges that define the perimeter of the region in the dual graph

        region_edges: Set[RegionTree.Edge] # Region tree edges connected to this vertex
        bridge_set: Set[TwinGraph.QuadEdge] # Edges that, if removed, would disconnect the primal tree
        bridge_to_region_edge_map: Dict[TwinGraph.QuadEdge, RegionTree.Edge] # Map from bridge edges to region tree edges

        point: Point

        # Instance labeling
        instance_counter: int = 0
        id_str: str

        def __init__(self, weight: float, dual_perimeter: List[Tuple[TwinGraph.QuadEdge, TwinGraph.EdgeDir]]) -> None:
            self.weight = weight
            self.dual_perimeter = dual_perimeter

            self.region_edges = set()
            self.bridge_set = set()
            self.bridge_to_region_edge_map = dict()

            self.id_str = f'R{RegionTree.Region.instance_counter:x}' # Hexadecimal ID for easier reading in debug
            RegionTree.Region.instance_counter += 1

            self.calc_point()

        def register_edge(self, edge: RegionTree.Edge) -> None:
            print("Registering edge", edge.twin_graph_edge.id_str, "to region", self.id_str)
            self.region_edges.add(edge)
            self.bridge_set.add(edge.twin_graph_edge)
            self.bridge_to_region_edge_map[edge.twin_graph_edge] = edge

        def unregister_edge(self, edge: RegionTree.Edge) -> None:
            print("Unregistering edge", edge.twin_graph_edge.id_str, "from region", self.id_str)
            self.region_edges.remove(edge)
            self.bridge_set.remove(edge.twin_graph_edge)
            del self.bridge_to_region_edge_map[edge.twin_graph_edge]

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
                # else:
                    # scaling -= 0.1 # Boost to push towards exterior

            if len(self.dual_perimeter) == 0:
                self.point = Point(0.0, 0.0)
            else:
                self.point = Point(x_sum / scaling, y_sum / scaling)

        # def get_all_dual_verts_contained_by_perimeter(self, graph: TwinGraph) -> Set[TwinGraph.Vert]:
        #     # Conduct a flood fill from one of the perimeter edges to find all dual verts to one
        #     #  side of the perimeter.  Then, check which side of the perimeter the external dual 
        #     #  vert is on and return the converse set of verts.
        #     if len(self.dual_perimeter) == 0:
        #         print("BBB")
        #         return graph.dualVerts

        #     # Get all verts in the perimeter
        #     perimeter_verts = self.get_perimeter_verts()
            
        #     # Use flood fill from the first perimeter edge's dual vert
        #     visited: Set[TwinGraph.Vert] = set()
        #     queue: List[TwinGraph.Vert] = []

        #     # Start from the dual vert of the first perimeter edge
        #     start_edge, _ = self.dual_perimeter[0]

        #     # Pick one of the dual verts (AB or BA)
        #     start_dual_vert = start_edge.dual_AB if start_edge.dual_AB is not None else start_edge.dual_BA
        #     if start_dual_vert is None:
        #         return set()
        #     queue.append(start_dual_vert)
        #     while queue:
        #         vert = queue.pop()
        #         if vert in visited:
        #             continue
        #         visited.add(vert)
                
        #         # Add all adjacent dual verts via cc_edges
        #         for edge in vert.cc_edges:
        #             # Get the other dual vert on this edge
        #             if edge.dual_AB == vert and edge.dual_BA is not None:
        #                 neighbor = edge.dual_BA
        #             elif edge.dual_BA == vert and edge.dual_AB is not None:
        #                 neighbor = edge.dual_AB
        #             else:
        #                 continue
        #             # Only add if not crossing the perimeter
        #             if neighbor not in perimeter_verts and neighbor not in visited:
        #                 queue.append(neighbor)

        #     print("!!!: ", graph.external_dual_vert in visited, len(visited))

        #     # If the external dual vert is in visited, we want the other side
        #     if graph.external_dual_vert in visited:
        #         # Return all dual verts not in visited
        #         return set(graph.dualVerts) - visited - perimeter_verts
        #     else:
        #         return visited - perimeter_verts
            
        def get_perimeter_verts(self) -> Set[TwinGraph.Vert]:
            verts = set()
            for edge, dir in self.dual_perimeter:
                src, _ = edge.get_dual_vert_pair(dir)
                verts.add(src)
            return verts

    class Edge:
        # Edge A and B ends must match directionally with underlying twin_graph_edge
        end_A: RegionTree.Region
        end_B: RegionTree.Region
        twin_graph_edge: TwinGraph.QuadEdge

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
