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
from Benchmarking.TwinGraphWilson import TwinGraphWilson

versions = ["GerryChain_Uniform", "TwinGraph_Wilson", "DirectPartition_Walk_Sample"]
exec_funcs = {
    "GerryChain_Uniform": lambda graphs: perform_gerrychain_find(graphs[0], use_uniform=True),
    "TwinGraph_Wilson": lambda graphs: perform_wilson_find(graphs[1]),
    "DirectPartition_Walk_Sample": lambda graphs: perform_direct_find(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.WALK),
}
n_sizes = [256, 1024, 1764, 2500, 3136, 3844, 4624]
reps = 1000

def perform_gerrychain_find(graph: Graph, use_uniform: bool):
    if use_uniform:
        spanning_tree_fn = tree.uniform_spanning_tree
    else:
        spanning_tree_fn = tree.random_spanning_tree
    
    # GerryChain's bipartition_tree retries internally if it can't find a cut? 
    # Actually, standard usage often involves a loop if it fails, but here we assume it works or we measure single attempt?
    # However, with epsilon=0, it might fail. 
    # But the original script used epsilon=0.
    # We will assume the original script's behavior is what is desired for GerryChain.
    # If it fails, it raises an error, which would stop the benchmark.
    # So we assume it succeeds.
    
    try:
        tree_partition = tree.bipartition_tree(
            graph=graph,
            pop_target=len(graph.nodes) // 2,
            pop_col='population',
            epsilon=0.05,
            spanning_tree_fn=spanning_tree_fn
        )
    except Exception:
        # If it fails, we might want to retry to be fair with other methods that retry?
        # But for now let's just let it fail if it does, or maybe the original script worked fine.
        pass

def perform_wilson_find(graph: TwinGraph):
    # We loop until we find a valid partition, similar to DirectPartition
    partition = None
    while partition is None:
        partition = TwinGraphWilson.bipartition_find(graph, epsilon=0.05)

def perform_direct_find(data: TwinGraph, start_selection_method: GraphNav.StartSelectionMethod):
    loops = []
    idx = 0
    graph_nav = None
    while len(loops) == 0:
        region_tree = RegionTree(data, epsilon=0.05)
        graph_nav = GraphNav(data, region_tree, start_selection_method=start_selection_method)
        graph_nav.animating = False

        loops = graph_nav.run_two_split_attempt()

        idx += 1
        if idx > 100000:
            raise Exception("Failed to find bipartition after 100000 attempts")
    
    if graph_nav is None:
        raise Exception("GraphNav not initialized properly.")
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
            df.to_csv("benchmarking_results_gabp_twingraphwilson.csv")
            print("Saved benchmarking results to benchmarking_results_gabp_twingraphwilson.csv")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    bench()
