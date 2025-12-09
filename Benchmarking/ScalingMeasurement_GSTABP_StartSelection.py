import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from math import sqrt, ceil
import time
import timeit

from wakepy import keep
import pandas as pd
import numpy as np
from gerrychain import Graph, tree

from TwinGraph import TwinGraph
from GraphNav import GraphNav
from RegionTree import RegionTree
from Euclid import *

from Benchmarking.GenerateGridGraphAjacencies import generate_grid_graph

versions = ["DirectPartition_Uniform", "DirectPartition_Central", "DirectPartition_Lining", "DirectPartition_Walk"]
exec_funcs = {
    "DirectPartition_Uniform": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.UNIFORM),
    "DirectPartition_Central": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.CENTRAL),
    "DirectPartition_Lining": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.LINING),
    "DirectPartition_Walk": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.WALK),
}
n_sizes = [256, 1024, 1764, 2500, 3136, 3844, 4624]
reps = 1000

def perform_direct_tree_attempt(data: TwinGraph, start_selection_method: GraphNav.StartSelectionMethod):
    region_tree = RegionTree(data, epsilon=0.05)
    graph_nav = GraphNav(data, region_tree, start_selection_method=start_selection_method)
    graph_nav.animating = False

    loops = graph_nav.run_two_split_attempt()
    if loops:
        partition_1 = graph_nav.get_enclosed_primal_verts(loops[0])
        partition_2 = set(data.primalVerts) - partition_1

def generate_n_grid_graphs(n: int) -> Tuple[Graph, TwinGraph]:
    # Validate input
    assert sqrt(n).is_integer(), "n must be a perfect square for grid graph generation."
    dim_n = int(sqrt(n))  # Ensure even dimensions

    # Generate Graph Source
    print(f"Generating TwinGraph for n={dim_n*dim_n} (dim={dim_n}x{dim_n})")
    points, weights, adjacencies = generate_grid_graph(dim_n, dim_n)

    # Create GerryChain graph
    gerrychain_graph = Graph(adjacencies)
    for n in gerrychain_graph.nodes:
        gerrychain_graph.nodes[n]['population'] = 1

    # Generate TwinGraph Source
    twin_graph = TwinGraph(points, weights, adjacencies)
    twin_graph.animating = False

    return gerrychain_graph, twin_graph

# @keep.running
def bench():
    with keep.running(): # Prevent sleep during benchmarking
        print("\n--- Starting Benchmarking ---\n")
        # Set up DF
        column_tuples = [(v, n) for v in versions for n in n_sizes]
        column_index = pd.MultiIndex.from_tuples(
            column_tuples,
            names=['version', 'n_size']
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

        for n in n_sizes:
            gerrychain_graph, twin_graph = generate_n_grid_graphs(n)
            graphs = (gerrychain_graph, twin_graph)

            for version in versions:
                print(f"Benchmarking {version} for n={n} over {reps} reps...")
                exec_func = exec_funcs[version]

                run_results = timeit.repeat(
                    stmt=lambda: exec_func(graphs),
                    timer=time.process_time,
                    repeat=reps,
                    number=1 # Acceptable to eat overhead since slowest runs > ~0.01s
                )
                df[(version, n)] = run_results

                avg_time = sum(run_results) / reps
                print(f"    Average Time: {avg_time:.6f} seconds")
        
            # Save results
            df.to_csv("benchmarking_results_gstabp_start_selection.csv")
            print("Saved benchmarking results to benchmarking_results_gstabp_start_selection.csv")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    bench()
