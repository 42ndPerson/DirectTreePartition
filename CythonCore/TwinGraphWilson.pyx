# cython: language_level=3

import random
from collections import deque
from CythonCore.TwinGraph cimport TwinGraph, Vert, QuadEdge

cdef class TwinGraphWilson:
    @staticmethod
    def bipartition_find(TwinGraph graph, double epsilon):
        """
        Runs Wilson's algorithm to build a uniform spanning tree on the TwinGraph,
        then searches for a cut that satisfies the epsilon balance constraint.
        Returns the two sets of vertices forming the partition, or None if no valid cut is found.
        """
        
        # 1. Build Uniform Spanning Tree using Wilson's Algorithm
        cdef dict tree_adj = TwinGraphWilson.wilson_spanning_tree(graph)
        
        # 2. Find a balanced cut
        cdef int total_weight = graph.total_primal_weight
        cdef double target_weight = total_weight / 2.0
        cdef double min_weight = target_weight * (1.0 - epsilon)
        cdef double max_weight = target_weight * (1.0 + epsilon)
        
        # Root the tree arbitrarily (e.g., at the first vertex) to compute subtree weights
        cdef Vert root = next(iter(graph.primalVerts))
        
        # Compute subtree weights and find valid cuts
        cdef dict subtree_weights = {}
        cdef dict parent_map = {root: None}
        
        # Post-order traversal stack
        cdef list stack = [root]
        cdef list dfs_order = []
        cdef Vert u, v
        
        while stack:
            u = stack.pop()
            dfs_order.append(u)
            for v in tree_adj[u]:
                if v != parent_map[u]:
                    parent_map[v] = u
                    stack.append(v)
        
        # Process in reverse DFS order (leaves first)
        cdef list valid_cuts = []
        cdef int weight
        
        for u in reversed(dfs_order):
            weight = u.weight
            for v in tree_adj[u]:
                if v != parent_map[u]: # v is a child of u
                    weight += subtree_weights[v]
            subtree_weights[u] = weight
            
            if u != root:
                if min_weight <= weight <= max_weight:
                    valid_cuts.append(u)
        
        if not valid_cuts:
            return None
            
        # 3. Return vertices partitioned by a chosen split edge
        cdef Vert cut_node = random.choice(valid_cuts)
        
        cdef set partition_a = set()
        cdef object q = deque([cut_node])
        
        while q:
            u = q.popleft()
            partition_a.add(u)
            for v in tree_adj[u]:
                if v != parent_map[u]: # Children
                    q.append(v)
                    
        cdef set partition_b = graph.primalVerts - partition_a
        
        return partition_a, partition_b

    @staticmethod
    def wilson_spanning_tree(TwinGraph graph):
        """
        Generates a uniform spanning tree using Wilson's algorithm.
        Returns an adjacency list representing the tree (undirected).
        """
        cdef list verts = list(graph.primalVerts)
        cdef Vert root = verts[0]
        cdef set in_tree = {root}
        
        cdef dict tree_adj = {v: [] for v in verts}
        
        random.shuffle(verts)
        
        cdef Vert u, curr, neighbor, next_v
        cdef QuadEdge edge
        cdef dict path
        
        for u in verts:
            if u in in_tree:
                continue
            
            curr = u
            path = {} # curr -> next
            
            while curr not in in_tree:
                # Pick random neighbor
                edge = random.choice(curr.cc_edges)
                neighbor = edge.get_primal_dest_from(curr)[0]
                
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
