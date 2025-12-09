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

versions = ["GerryChain_Uniform", "GerryChain_MST", "DirectPartition_Walk_Sample"]
exec_funcs = {
    "GerryChain_Uniform": lambda graphs: perform_gerrychain_find(graphs[0], use_uniform=True),
    "GerryChain_MST": lambda graphs: perform_gerrychain_find(graphs[0], use_uniform=False),
    "DirectPartition_Walk_Sample": lambda graphs: perform_direct_find(graphs[1], start_selection_method=GraphNav.StartSelectionMethod.WALK),
}
# Shapes: (width, height)
# We aim for roughly 2500 nodes, with increasing aspect ratios.
# Since epsilon is 0, we ensure even number of nodes.
shapes = [
    (50, 50),   # 2500, 1:1
    (36, 70),   # 2520, 1:2
    (25, 100),  # 2500, 1:4
    (18, 140),  # 2520, 1:8
    (12, 208),  # 2496, 1:16
    (9, 280),   # 2520, 1:31
    (6, 416),   # 2496, 1:69
    (4, 625)    # 2500, 1:156
]
# Pre-process shapes to ensure even node count
adjusted_shapes = []
for w, h in shapes:
    if (w * h) % 2 != 0:
        h += 1
    adjusted_shapes.append((w, h))
shapes = adjusted_shapes

reps = 1000

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
    loops = []
    idx = 0
    graph_nav = None
    while not loops:
        region_tree = RegionTree(data)
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

def generate_rect_grid_graphs(w: int, h: int) -> Tuple[Graph, TwinGraph]:
    # Generate Graph Source
    print(f"Generating TwinGraph for {w}x{h} (n={w*h})")
    points, weights, adjacencies = generate_grid_graph(w, h)

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
            df.to_csv("benchmarking_results_rgebp.csv")
            print("Saved benchmarking results to benchmarking_results_rgebp.csv")

        print("\n--- Benchmarking Results ---")
        print(df.describe().T)

if __name__ == "__main__":
    # Note: It is essential to run this script as user initiated.
    #  This should lead to a high Quality of Service (QoS) on macOS,
    #  preventing the system from sending to E-cores.
    bench()