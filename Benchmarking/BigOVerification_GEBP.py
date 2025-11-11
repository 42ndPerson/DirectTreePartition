import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from math import sqrt, ceil

import big_o

from TwinGraph import TwinGraph
from GraphNav import GraphNav
from RegionTree import RegionTree
from Euclid import *

from Benchmarking.GenerateGridGraphAjacencies import generate_grid_graph

# This code is meant to validate the Big O performance of the bipartition finding algorithm
# Note: big_o library uses wall clock time, so results may vary based on system load

def perform_find(data: TwinGraph):
    loop = None
    idx = 0
    graph_nav = None
    while loop == None:
        region_tree = RegionTree(data)
        graph_nav = GraphNav(data, region_tree)
        graph_nav.animating = False

        loop = graph_nav.run_two_split_attempt()

        idx += 1
        # print("        Attempt:", idx)
        if idx > 100000:
            raise Exception("Failed to find bipartition after 100000 attempts")
    # print(f"    Found bipartition after {idx} attempts")
    
    if graph_nav is None:
        raise Exception("GraphNav not initialized properly.")
    partition_1 = graph_nav.get_enclosed_primal_verts(loop)
    partition_2 = set(data.primalVerts) - partition_1
    return True

def generate_n_grid_twingraph(n: int) -> TwinGraph:
    dim_n = int(ceil(sqrt(n) / 2) * 2)  # Ensure even dimensions
    print(f"Generating TwinGraph for n={dim_n*dim_n} (dim={dim_n}x{dim_n})")
    points, weights, adjacencies = generate_grid_graph(dim_n, dim_n)
    return TwinGraph(points, weights, adjacencies)

best_fit, all_fits = big_o.big_o(perform_find, generate_n_grid_twingraph, n_repeats=30, max_n=16384)

print("--- Best Fit ---")
print(best_fit)

print("\n--- All Fits ---")
for class_name, residual in all_fits.items():
    print(f"{class_name}: {residual:.2e}")

