from typing import List, Tuple

from CythonCore.Euclid import Point
from CythonCore.TwinGraph import TwinGraph
from GraphVisCython import GraphVis
from CythonCore.RegionTree import RegionTree
from CythonCore.GraphNav import GraphNav

def generate_grid_graph(width: int, height: int) -> Tuple[List[Point], List[int], List[Tuple[int, int]]]:
    points = []
    weights = []
    edgeIdxs = []

    for row in range(height):
        for column in range(width):
            idx = row*width + column
            points.append(Point(column, row))
            weights.append(1)
            if column != width-1:
                edgeIdxs.append((idx, idx+1))
            if row != height-1:
                edgeIdxs.append((idx, idx+width))

    return points, weights, edgeIdxs

points, weights, edge_idxs = generate_grid_graph(70,70)
graph = TwinGraph(points, weights, edge_idxs)
region_tree = RegionTree(graph, epsilon=0)
graph_nav = GraphNav(graph, region_tree, multi_walk_attempts=15)

GraphVis(graph, graph_nav, region_tree, 0)
