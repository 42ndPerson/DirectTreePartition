from typing import List, Tuple

from Euclid import *
from TwinGraph import *
from GraphVis import *
from RegionTree import *

def generate_grid_graph(width: int, height: int) -> Tuple[List[Point], List[float], List[Tuple[int, int]]]:
    points = []
    weights = []
    edgeIdxs = []

    for row in range(height):
        for column in range(width):
            idx = row*width + column
            points.append(Point(column, row))
            weights.append(1)#idx)
            if column != width-1:
                edgeIdxs.append((idx, idx+1))
            if row != height-1:
                edgeIdxs.append((idx, idx+width))

    return points, weights, edgeIdxs

points, weights, edge_idxs = generate_grid_graph(26,26)
graph = TwinGraph(points, weights, edge_idxs)
region_tree = RegionTree(graph)
graph_nav = GraphNav(graph, region_tree)

GraphVis(graph, graph_nav, region_tree)