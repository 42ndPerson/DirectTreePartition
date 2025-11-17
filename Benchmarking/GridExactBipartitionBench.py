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
from Benchmarking.GenerateGridGraphAjacencies import generate_grid_graph


def benchmark_bipartitioning_for_g_grid_graph(g: int, iters: int):
    points, weights, adjacencies = generate_grid_graph(g, g)

    # Create GerryChain graph
    gerrychain_graph = Graph(adjacencies)
    for n in gerrychain_graph.nodes:
        gerrychain_graph.nodes[n]['population'] = 1
    gerrychain_target_population = len(gerrychain_graph.nodes) / 2

    # Create direct partition graph
    graph = TwinGraph(points, weights, adjacencies)
    graph.animating = False

    # Run Gerrychain
    # start = time.perf_counter()
    # for i in range(iters):
    #     tree_partition = tree.bipartition_tree(
    #         graph=gerrychain_graph,
    #         pop_target=gerrychain_target_population,
    #         pop_col='population',
    #         epsilon=0,
    #         spanning_tree_fn=tree.uniform_spanning_tree
    #         # spanning_tree_fn=tree.random_spanning_tree
    #     )
    #     print("GerryChain Tree:", i)
    # end = time.perf_counter()
    # print(f"GerryChain Bipartitioning Time for {iters} iterations: {end - start:.4f} seconds")
    
    # Run direct partitioning
    start = time.perf_counter()
    for i in range(iters):
        loop = None
        idx = 0
        graph_nav = None
        while loop == None:
            region_tree = RegionTree(graph)
            graph_nav = GraphNav(graph, region_tree)
            graph_nav.animating = False

            loop = graph_nav.run_two_split_attempt()

            idx += 1
            if idx > 100000:
                raise Exception("Failed to find bipartition after 100000 attempts")
            
        if graph_nav is None:
            raise Exception("GraphNav not initialized properly.")
        partition_1 = graph_nav.get_enclosed_primal_verts(loop)
        partition_2 = set(graph.primalVerts) - partition_1
        print("Direct Bipartition:", i, "after", idx, "attempts")
    end = time.perf_counter()
    print(f"Direct Bipartitioning Time for {iters} iterations: {end - start:.4f} seconds")

profiler = cProfile.Profile()
with profiler:  
    benchmark_bipartitioning_for_g_grid_graph(100, 100)

profiler.disable()
profiler.dump_stats('profile.stats')