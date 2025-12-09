import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from math import sqrt, ceil
import time
import timeit
from typing import List, Tuple, Set, Optional

from wakepy import keep
import pandas as pd
import numpy as np
from gerrychain import Graph, tree

from CythonCore.TwinGraph import TwinGraph
from CythonCore.GraphNav import GraphNav, StartSelectionMethod
from CythonCore.RegionTree import RegionTree
from CythonCore.Euclid import Point
from CythonCore.TwinGraphWilson import TwinGraphWilson

from Benchmarking.GenerateGridGraphAjacencies import generate_grid_graph

versions = ["DirectPartition_Cython", "TwinGraph_Wilson_Cython"]
exec_funcs = {
    "DirectPartition_Cython": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=StartSelectionMethod.WALK),
    "TwinGraph_Wilson_Cython": lambda graphs: perform_wilson_tree_attempt(graphs[1]),
}

# Define target sizes and aspect ratios
base_ns = [2500, 4900, 8100]
ratios = [1, 2, 4, 8, 16, 32]

shapes = []
for n in base_ns:
    for r in ratios:
        # h / w = r  => h = r*w
        # w * h = n  => w * r*w = n => w^2 = n/r => w = sqrt(n/r)
        w = int(sqrt(n / r))
        h = int(n / w)
        # Ensure even node count for bipartition
        if (w * h) % 2 != 0:
            h += 1
        shapes.append((w, h))

reps = 1000

def perform_wilson_tree_attempt(graph: TwinGraph):
    TwinGraphWilson.wilson_spanning_tree(graph)

def perform_direct_tree_attempt(data: TwinGraph, start_selection_method: StartSelectionMethod):
    region_tree = RegionTree(data)
    graph_nav = GraphNav(data, region_tree, start_selection_method=start_selection_method)
    graph_nav.animating = False

    loops = graph_nav.run_two_split_attempt()
    if loops:
        partition_1 = graph_nav.get_enclosed_primal_verts(loops[0])
        partition_2 = set(data.primalVerts) - partition_1

def generate_rect_grid_graphs(w: int, h: int) -> Tuple[Graph, TwinGraph]:
    # Generate Graph Source
    print(f"Generating TwinGraph for {w}x{h} (n={w*h})")
    points, weights, adjacencies = generate_grid_graph(w, h)

    # Create GerryChain graph
    gerrychain_graph = Graph(adjacencies)
    for n in gerrychain_graph.nodes:
        gerrychain_graph.nodes[n]['population'] = 1

    # Generate TwinGraph Source
    # Convert points to CythonCore.Euclid.Point objects
    cython_points = [Point(p.x, p.y) for p in points]
    twin_graph = TwinGraph(cython_points, weights, adjacencies)
    twin_graph.animating = False

    return gerrychain_graph, twin_graph

# @keep.running
def bench():
    with keep.running(): # Prevent sleep during benchmarking
        print("\n--- Starting Benchmarking (Cython Rectangular) ---\n")
        # Set up DF
        shape_strs = [f"{w}x{h}" for w, h in shapes]
        column_tuples = [(v, s) for v in versions for s in shape_strs]
        column_index = pd.MultiIndex.from_tuples(
            column_tuples,
            names=['version', 'shape']
        )

        # Index run number
        row_index = pd.Index(range(reps), name='run')

        # Create the empty DataFrame
        # Pre-fill with np.nan
        df = pd.DataFrame(
            np.nan,
            index=row_index,
            columns=column_index
        )

        for w, h in shapes:
            shape_str = f"{w}x{h}"
            gerrychain_graph, twin_graph = generate_rect_grid_graphs(w, h)
            graphs = (gerrychain_graph, twin_graph)

            for version in versions:
                print(f"Benchmarking {version} for {shape_str} over {reps} reps...")
                exec_func = exec_funcs[version]

                run_results = timeit.repeat(
                    stmt=lambda: exec_func(graphs),
                    timer=time.process_time,
                    repeat=reps,
                    number=1 # Acceptable to eat overhead since slowest runs > ~0.01s
                )
                df[(version, shape_str)] = run_results

                avg_time = sum(run_results) / reps
                print(f"    Average Time: {avg_time:.6f} seconds")
        
            # Save results
            df.to_csv("benchmarking_results_rgstebp_cython.csv")
            print("Saved benchmarking results to benchmarking_results_rgstebp_cython.csv")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    # Note: It is essential to run this script as user initiated.
    #  This should lead to a high Quality of Service (QoS) on macOS,
    #  preventing the system from sending to E-cores.
    bench()
