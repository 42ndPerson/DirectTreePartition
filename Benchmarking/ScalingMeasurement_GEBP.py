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

versions = ["GerryChain_Uniform", "GerryChain_MST", "DirectPartition_Uniform", "DirectPartition_Central"]
exec_funcs = {
    "GerryChain_Uniform": lambda graphs: perform_gerrychain_find(graphs[0], use_uniform=True),
    "GerryChain_MST": lambda graphs: perform_gerrychain_find(graphs[0], use_uniform=False),
    "DirectPartition_Uniform": lambda graphs: perform_direct_find(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.UNIFORM),
    "DirectPartition_Central": lambda graphs: perform_direct_find(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.CENTRAL),
}
n_sizes = [256, 1024, 4096]#, 29241, 58564, 87616, 116964, 145924, 174724, 204304, 233289, 262144]
reps = 100

def perform_gerrychain_find(graph: Graph, use_uniform: bool):
    if use_uniform:
        spanning_tree_fn = tree.uniform_spanning_tree
    else:
        spanning_tree_fn = tree.random_spanning_tree
    tree_partition = tree.bipartition_tree(
        graph=graph,
        pop_target=len(graph.nodes) // 2,
        pop_col='population',
        epsilon=0,
        spanning_tree_fn=spanning_tree_fn
    )

def perform_direct_find(data: TwinGraph, start_selection_method: GraphNav.StartSelectionMethod):
    loop = None
    idx = 0
    graph_nav = None
    while loop == None:
        region_tree = RegionTree(data)
        graph_nav = GraphNav(data, region_tree, start_selection_method=start_selection_method)
        graph_nav.animating = False

        loop = graph_nav.run_two_split_attempt()

        idx += 1
        if idx > 100000:
            raise Exception("Failed to find bipartition after 100000 attempts")
    
    if graph_nav is None:
        raise Exception("GraphNav not initialized properly.")
    partition_1 = graph_nav.get_enclosed_primal_verts(loop)
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

    return gerrychain_graph, twin_graph

# @keep.running
def bench():
    with keep.running():
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
                    repeat=reps,
                    number=1 # Acceptable to eat overhead since slowest runs > ~0.01s
                )
                df[(version, n)] = run_results

                avg_time = sum(run_results) / reps
                print(f"    Average Time: {avg_time:.6f} seconds")
        
        # Save results
        df.to_csv("benchmarking_results.csv")
        print("Saved benchmarking results to benchmarking_results.csv")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    bench()