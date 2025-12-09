import sys
import os
import random
from collections import deque
from typing import Set, Tuple, Dict, List, Optional

# Adjust path to allow imports from parent directory
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from TwinGraph import TwinGraph

class TwinGraphWilson:
    @staticmethod
    def bipartition_find(graph: TwinGraph, epsilon: float) -> Optional[Tuple[Set[TwinGraph.Vert], Set[TwinGraph.Vert]]]:
        """
        Runs Wilson's algorithm to build a uniform spanning tree on the TwinGraph,
        then searches for a cut that satisfies the epsilon balance constraint.
        Returns the two sets of vertices forming the partition, or None if no valid cut is found.
        """
        
        # 1. Build Uniform Spanning Tree using Wilson's Algorithm
        tree_adj = TwinGraphWilson.wilson_spanning_tree(graph)
        
        # 2. Find a balanced cut
        total_weight = graph.total_primal_weight
        target_weight = total_weight / 2.0
        min_weight = target_weight * (1 - epsilon)
        max_weight = target_weight * (1 + epsilon)
        
        # Root the tree arbitrarily (e.g., at the first vertex) to compute subtree weights
        root = next(iter(graph.primalVerts))
        
        # Compute subtree weights and find valid cuts
        # We need to traverse the tree.
        # parent_map: child -> parent (to avoid traversing back up)
        # subtree_weights: vert -> weight
        
        subtree_weights: Dict[TwinGraph.Vert, int] = {}
        parent_map: Dict[TwinGraph.Vert, Optional[TwinGraph.Vert]] = {root: None}
        
        # Post-order traversal stack
        # We can do a topological sort or just a DFS to build the stack
        stack = [root]
        dfs_order = []
        while stack:
            u = stack.pop()
            dfs_order.append(u)
            for v in tree_adj[u]:
                if v != parent_map[u]:
                    parent_map[v] = u
                    stack.append(v)
        
        # Process in reverse DFS order (leaves first)
        valid_cuts = []
        
        for u in reversed(dfs_order):
            weight = u.weight
            for v in tree_adj[u]:
                if v != parent_map[u]: # v is a child of u
                    weight += subtree_weights[v]
            subtree_weights[u] = weight
            
            # Check if cutting the edge above u (between u and parent_map[u]) is valid
            # The component defined by u has weight 'weight'.
            # The other component has 'total_weight - weight'.
            # We only need to check 'weight' against bounds because total_weight is fixed.
            
            if u != root: # Root has no edge above it
                if min_weight <= weight <= max_weight:
                    valid_cuts.append(u)
        
        if not valid_cuts:
            return None
            
        # 3. Return vertices partitioned by a chosen split edge
        # If multiple cuts are valid, GerryChain typically chooses uniformly at random?
        # Or just picks one. We'll pick one randomly.
        cut_node = random.choice(valid_cuts)
        
        # Reconstruct the partition sets
        # Partition A: The subtree rooted at cut_node (in our directed view)
        # Partition B: The rest
        
        partition_a = set()
        q = deque([cut_node])
        while q:
            u = q.popleft()
            partition_a.add(u)
            for v in tree_adj[u]:
                if v != parent_map[u]: # Children
                    q.append(v)
                    
        partition_b = graph.primalVerts - partition_a
        
        return partition_a, partition_b

    @staticmethod
    def wilson_spanning_tree(graph: TwinGraph) -> Dict[TwinGraph.Vert, List[TwinGraph.Vert]]:
        """
        Generates a uniform spanning tree using Wilson's algorithm.
        Returns an adjacency list representing the tree (undirected).
        """
        verts = list(graph.primalVerts)
        root = verts[0] # Arbitrary root for the tree construction
        in_tree = {root}
        
        # Adjacency list for the spanning tree
        tree_adj: Dict[TwinGraph.Vert, List[TwinGraph.Vert]] = {v: [] for v in verts}
        
        # Shuffle verts to process in random order (optimization for Wilson's)
        # We iterate through all vertices to ensure they are connected to the tree
        random.shuffle(verts)
        
        for u in verts:
            if u in in_tree:
                continue
            
            # Perform loop-erased random walk from u until we hit the tree
            curr = u
            path: Dict[TwinGraph.Vert, TwinGraph.Vert] = {} # curr -> next
            
            while curr not in in_tree:
                # Pick random neighbor
                # TwinGraph.Vert.cc_edges contains QuadEdges
                # We need to get the neighbor vertex
                edge = random.choice(curr.cc_edges)
                neighbor, _ = edge.get_primal_dest_from(curr)
                
                path[curr] = neighbor
                curr = neighbor
            
            # Add the loop-erased path to the tree
            curr = u
            while curr not in in_tree:
                next_v = path[curr]
                in_tree.add(curr)
                tree_adj[curr].append(next_v)
                tree_adj[next_v].append(curr)
                curr = next_v
                
        return tree_adj
