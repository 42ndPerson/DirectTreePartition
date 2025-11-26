import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from math import sqrt
import timeit
import time

from wakepy import keep
import pandas as pd
import numpy as np
from gerrychain import Graph, tree

from TwinGraph import TwinGraph
from GraphNav import GraphNav
from RegionTree import RegionTree
from Euclid import *

from Benchmarking.GenerateGridGraphAjacencies import generate_grid_graph
from typing import Tuple

# Configuration
versions = ["GerryChain_Uniform", "GerryChain_MST", "DirectPartition_Walk_Single"]
n_sizes = [256, 1296, 2500, 3600, 4624, 5776, 6724, 7744, 8836, 10000, 250000] 
reps = 100

def perform_gerrychain_tree_attempt(graph: Graph, use_uniform: bool):
    """
    Generates a single spanning tree AND evaluates it for valid cuts.
    """
    if use_uniform:
        spanning_tree_fn = tree.uniform_spanning_tree
    else:
        spanning_tree_fn = tree.random_spanning_tree

    try:
        tree_partition = tree.bipartition_tree(
            graph=graph,
            pop_target=len(graph.nodes) // 2,
            pop_col='population',
            epsilon=0,
            spanning_tree_fn=spanning_tree_fn,
            max_attempts=1  # Only attempt once
        )
    except Exception:
        pass  # Ignore failures for benchmarking

def perform_direct_tree_attempt(data: TwinGraph, start_selection_method: GraphNav.StartSelectionMethod):
    """
    Performs exactly one attempt at building a RegionTree to find a split.
    Does NOT retry if the attempt fails (returns None).
    """
    region_tree = RegionTree(data)
    graph_nav = GraphNav(data, region_tree, start_selection_method=start_selection_method)
    graph_nav.animating = False

    # run_two_split_attempt performs the walks and region divisions.
    # It inherently calculates weights and checks for splits (validate_split_possible).
    # It returns a loop if successful, or None if it hits a dead end.
    graph_nav.run_two_split_attempt()

# Execution mappings
exec_funcs = {
    "GerryChain_Uniform": lambda graphs: perform_gerrychain_tree_attempt(graphs[0], use_uniform=True),
    "GerryChain_MST": lambda graphs: perform_gerrychain_tree_attempt(graphs[0], use_uniform=False),
    "DirectPartition_Walk_Single": lambda graphs: perform_direct_tree_attempt(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.WALK),
}

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
        print("\n--- Starting Single-Tree Exploration Benchmarking ---\n")
        # Set up DF
        column_tuples = [(v, n) for v in versions for n in n_sizes]
        column_index = pd.MultiIndex.from_tuples(
            column_tuples,
            names=['version', 'n_size']
        )

        # Index run number
        row_index = pd.Index(range(reps), name='run')

        # Create the empty DataFrame
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
                    number=1 
                )
                df[(version, n)] = run_results

                avg_time = sum(run_results) / reps
                print(f"    Average Time: {avg_time:.6f} seconds")
        
            # Save results
            output_file = "benchmarking_tree_explore_results.csv"
            df.to_csv(output_file)
            print(f"Saved benchmarking results to {output_file}")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    bench()