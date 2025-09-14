from Euclid import *
from HyperGraph import *
from GraphVis import *

def generate_grid_graph(width: int, height: int) -> ([Point], [float], [(int, int)]):
    points = []
    weights = []
    edgeIdxs = []

    for row in range(height):
        for column in range(width):
            idx = row*width + column
            points.append(Point(column, row))
            weights.append(idx)
            if column != width-1:
                edgeIdxs.append((idx, idx+1))
            if row != height-1:
                edgeIdxs.append((idx, idx+width))

    return points, weights, edgeIdxs

points, weights, edge_idxs = generate_grid_graph(25,25)
graph = HyperGraph(points, weights, edge_idxs)

GraphVis(graph)