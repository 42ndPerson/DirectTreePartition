from typing import List, Tuple

from Euclid import *
from TwinGraph import *
from GraphVis import *
from RegionTree import *
import networkx as nx
import matplotlib.pyplot as plt

# Quick way of drawing graphs using wasd as arrow keys to move around a grid and space to draw an edge.
def wasd_code_graph(s):
    x = 0
    y = 0
    last = (x, y)
    m = {last: 0}
    next_id = 1
    jump = False
    edges = set()
    for c in s:
        if c == "a":
            x -= 1
        elif c == "d":
            x += 1
        elif c == "s":
            y -= 1
        elif c == "w":
            y += 1
        elif c == "j":
            jump = True
        elif c == " ":
            if (x, y) not in m:
                m[(x, y)] = next_id
                next_id += 1
            if jump:
                jump = False
            else:
                edges.add((m[last], m[(x, y)]))
            last = (x, y)
    return edges


# Lats and Lons from https://www.ala.org/magirt/publicationsab/ia
def generate_iowa_graph() -> Tuple[List[Point], List[int], List[Tuple[int, int]]]:
    points = [None for _ in range(99)]
    weights = [None for _ in range(99)]

    with open("Data/IACountiesLatLon.txt", "r") as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        a = line.split(" ")
        if len(a) != 11:
            raise Exception()
        idx = int(a[9])
        column = (int(a[5]) + int(a[6])/60 + int(a[7]) + int(a[8])/60)/2
        row = (int(a[1]) + int(a[2])/60 + int(a[3]) + int(a[4])/60)/2
        points[idx] = Point(column, row)
        weights[idx] = int(a[10])

    edgeIdxs = wasd_code_graph("s s s s s jwwwwwd a d s a d s a d sa dw s a d s a d sa dw s a d s aw sd s s jwwwwwwwwd a d s a d s a d s a d sa dw s a d s a d s aw sd a d s wa ds a d s a d jwwwwwwwwd a d s a d s a d s a d sa dw s a d s aw sd a d s a d s aw sd a d s a d jwwwwwwwd wa s d s a d s aw sd a d sa dw s a d s aw sd a d s a d s aw sd a d s a d jwwwwwwwwd as d w s s a d sa dw s a d sa dw s a d s aw sd a d s a d s aw sd a d s a d jwwwwwwwwd a d s a d s a d s a d sa dw s a d s aw sd a d s aw sd a d s aw sd a d s a d jwwwwwwwwd a s dw s a d s a d s a d sa dw s a d s a d s aw sd a d s aw sd a d s a d jwwwwwwwwd a d s a d s a d s aw s d sa dw s a d s a d s aw sd a d s aw sd a d s a d jwwwwwwwd wa s d s aw s d s a d sa dw s a d s a d s aw sd a d s aw sd a d s a d jwwwwwwwd a d s a d s a d sa dw s a d s aw sd a d s wa ds s wwa dss wa ds s wwa dss wa ds a d jwwwwwd wa s d sa dw s a d s wa ds a d s wa ds a s d sa dw s a d aa ")
    return points, weights, list(edgeIdxs)

points, weights, edge_idxs = generate_iowa_graph()

if False:
    # Visualize graph in networkx
    G = nx.Graph()
    for i, point in enumerate(points):
        G.add_node(i, pos=(-point.y, point.x))

    for u, v in edge_idxs:
        G.add_edge(u, v)

    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, with_labels=True, node_size=50, font_size=8)
    plt.savefig('Plots/iowa_graph.png', dpi=150, bbox_inches='tight')
    print("Graph saved to Plots/iowa_graph.png")
else:
    graph = TwinGraph(points, weights, edge_idxs)
    region_tree = RegionTree(graph)
    graph_nav = GraphNav(graph, region_tree)

    GraphVis(graph, graph_nav, region_tree)