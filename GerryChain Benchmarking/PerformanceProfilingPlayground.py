# Import modules from enclosing driectory
import sys
import os
import time
import cProfile
from typing import List, Tuple
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from gerrychain import Graph, tree
from matplotlib import pyplot as plt

from TwinGraph import TwinGraph
from GraphNav import GraphNav
from RegionTree import RegionTree
from Euclid import *

import numpy as np
def generate_grid_graph(width: int, height: int) -> Tuple[List[Point], List[int], List[Tuple[int, int]]]:
    w, h = int(width), int(height)
    n = w * h

    # Points
    xs = np.tile(np.arange(w, dtype=int), h)
    ys = np.repeat(np.arange(h, dtype=int), w)
    points = [Point(int(x), int(y)) for x, y in zip(xs, ys)]

    weights = [1] * n

    # Horizontal edges
    if w > 1 and h > 0:
        start_h = (np.arange(h, dtype=int)[:, None] * w + np.arange(w - 1, dtype=int)[None, :]).ravel()
        end_h = start_h + 1
        horiz = np.stack((start_h, end_h), axis=1)
    else:
        horiz = np.empty((0, 2), dtype=int)

    # Vertical edges
    if h > 1 and w > 0:
        start_v = (np.arange(h - 1, dtype=int)[:, None] * w + np.arange(w, dtype=int)[None, :]).ravel()
        end_v = start_v + w
        vert = np.stack((start_v, end_v), axis=1)
    else:
        vert = np.empty((0, 2), dtype=int)

    edges = np.vstack((horiz, vert))
    edgeIdxs = [(int(a), int(b)) for a, b in edges]

    return points, weights, edgeIdxs

def benchmark_for_g_grid_graph(g: int):
    points, weights, adjecencies = generate_grid_graph(g, g)

    # Create GerryChain graph
    gerrychain_graph = Graph(adjecencies)
    for n in gerrychain_graph.nodes:
        gerrychain_graph.nodes[n]['population'] = 1
    gerrychain_target_population = len(gerrychain_graph.nodes) / 2

    # Create direct partition graph
    graph = TwinGraph(points, weights, adjecencies)
    graph.animating = False

    profiler = cProfile.Profile()
    with profiler:
        start = time.perf_counter()
        edge: TwinGraph.QuadEdge = next(iter(graph.edges))
        vert, _ = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
        for i in range(100000):
            res = edge.get_dual_cc_next_edge(vert)
        end = time.perf_counter()
        print(f"QuadEdge.get_dual_cc_next_edge Time for 10000 iterations: {end - start:.4f} seconds")

        start = time.perf_counter()
        edge: TwinGraph.QuadEdge = next(iter(graph.edges))
        vert, _ = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
        for i in range(100000):
            res = None
            if vert == edge.dual_AB:
                _, dir = None, TwinGraph.EdgeDir.AB #edge.dual_AB_cc_next.get_dual_dest_from(vert)
                res = (edge.dual_AB_cc_next, dir)
            if vert == edge.dual_BA:
                _, dir = None, TwinGraph.EdgeDir.BA #edge.dual_BA_cc_next.get_dual_dest_from(vert)
                res = (edge.dual_BA_cc_next, dir)
        end = time.perf_counter()
        print(f"QuadEdge direct conditional Time for 10000 iterations: {end - start:.4f} seconds")

        start = time.perf_counter()
        edge: TwinGraph.QuadEdge = next(iter(graph.edges))
        vert, _ = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
        for i in range(100000):
            res = edge.dual_cc_next_retrieval_cache[id(vert)]
        end = time.perf_counter()
        print(f"QuadEdge cc lookup cache access Time for 10000 iterations: {end - start:.4f} seconds")

        start = time.perf_counter()
        edge: TwinGraph.QuadEdge = next(iter(graph.edges))
        vert, _ = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
        for i in range(100000):
            res = edge.dual_AB_cc_next
        end = time.perf_counter()
        print(f"QuadEdge cc direct access Time for 10000 iterations: {end - start:.4f} seconds")

    profiler.disable()
    profiler.dump_stats('profile.stats')

benchmark_for_g_grid_graph(100)