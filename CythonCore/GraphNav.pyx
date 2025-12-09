# cython: language_level=3
import random
from collections import deque
from CythonCore.TwinGraph cimport TwinGraph, Vert, QuadEdge, EdgeDir, VertRole, reverse_dir
from CythonCore.RegionTree cimport RegionTree, Region, Edge

cdef class GraphNav:
    def __init__(self, TwinGraph graph, RegionTree region_tree, 
                 StartSelectionMethod start_selection_method=StartSelectionMethod.WALK,
                 MultiWalkStartBehavior multi_walk_start_behavior=MultiWalkStartBehavior.HIT_POINT,
                 int multi_walk_attempts=5,
                 bint unwind_loops=False):
        self.graph = graph
        self.region_tree = region_tree
        self.target_regions = deque(region_tree.regions)
        
        self.tree_verts = set()
        self.tree_edges = set()
        self.bridge_edges = set()
        
        self.start_selection_method = start_selection_method
        self.multi_walk_start_behavior = multi_walk_start_behavior
        self.multi_walk_attempts = multi_walk_attempts
        self.unwind_loops = unwind_loops
        
        self.animating = True
        self.animation_tracks = [[[]], [[]]]

    cpdef list run_two_split_attempt(self):
        cdef Region target_region
        cdef Vert origin_vert = None
        cdef set lining_verts
        cdef set new_target_regions
        cdef list loops = []
        cdef list loop
        cdef Edge edge_center
        cdef list loop_edges
        cdef set seen
        cdef set duplicates
        cdef list boundary
        cdef QuadEdge edge
        
        while self.target_regions:
            target_region = self.target_regions.popleft()
            
            if self.start_selection_method == StartSelectionMethod.UNIFORM:
                origin_vert = target_region.get_uniform_interior_vert()
            elif self.start_selection_method == StartSelectionMethod.CENTRAL:
                origin_vert = target_region.get_central_region_vert()
            elif self.start_selection_method == StartSelectionMethod.LINING:
                lining_verts, _ = target_region.get_interior_lining_verts()
                if len(lining_verts) == 0:
                    continue
                origin_vert = random.choice(list(lining_verts))
            elif self.start_selection_method == StartSelectionMethod.WALK:
                origin_vert = target_region.get_walk_random_interior_vert(steps=(len(target_region.dual_perimeter)//1))
                
            if origin_vert is None:
                continue
                
            new_target_regions = self.walk_division_from(target_region, origin_vert)
            self.target_regions.extend(new_target_regions)
            
        if len(self.region_tree.edge_centers) > 0:
            for edge_center in self.region_tree.edge_centers:
                loop = self.traverse_clockwise_loop(edge_center.twin_graph_edge, EdgeDir.AB, {edge_center.twin_graph_edge})
                
                if self.animating:
                    loop_edges = [e for e, _ in loop]
                    seen = set()
                    duplicates = {edge_center.twin_graph_edge}
                    for edge in loop_edges:
                        if edge in seen:
                            duplicates.add(edge)
                        else:
                            seen.add(edge)
                    loop_edges = [e for e in loop_edges if e not in duplicates]
                    
                    boundary = [
                        (edge, VertRole.DUAL, EdgeDir.AB, 0)
                        for edge in loop_edges
                    ]
                    boundary.append((edge_center.twin_graph_edge, VertRole.DUAL, EdgeDir.AB, 1))
                    self.animation_tracks[0].append(boundary)
                    
                loops.append(loop)
            return loops
        else:
            return []

    cpdef set walk_division_from(self, Region region_in, Vert origin_vert):
        cdef Region region = None
        cdef set perimeter_verts
        cdef Region candidate_region
        
        if region_in is None:
            perimeter_verts = self.get_enclosing_perimeter_verts(origin_vert)
            for candidate_region in self.region_tree.regions:
                if candidate_region.get_perimeter_verts().issuperset(perimeter_verts):
                    region = candidate_region
                    break
            if region is None:
                raise ValueError("Starting vert is not contained in any region.")
        else:
            region = region_in
            
        cdef list walk_edges = []
        cdef set consumed_verts = set()
        cdef set new_target_regions = set()
        cdef Vert start_vert = origin_vert
        cdef int i
        cdef list new_walk_edges
        cdef set new_consumed_verts
        cdef Vert termination_vert
        cdef QuadEdge region_edge
        cdef list new_regions
        cdef Region new_region
        
        for i in range(self.multi_walk_attempts):
            new_walk_edges, new_consumed_verts, termination_vert = self.loop_erased_random_walk_from(start_vert, walk_edges, consumed_verts, region)
            if self.multi_walk_start_behavior == MultiWalkStartBehavior.HIT_POINT:
                if termination_vert not in self.tree_verts and termination_vert.role != VertRole.DUAL_EXTERIOR:
                    start_vert = termination_vert
                else:
                    start_vert = origin_vert
            
            walk_edges.extend(new_walk_edges)
            consumed_verts.update(new_consumed_verts)
            
        region_edge = self.commit_walk(walk_edges, consumed_verts)
        new_regions = self.develop_region(region, region_edge, EdgeDir.AB)
        
        for new_region in new_regions:
            new_region.generate_cc_region_edges()
        for new_region in new_regions:
            if new_region.check_center():
                new_target_regions.add(new_region)
                
        self.region_tree.remove_region(region)
        return new_target_regions

    cpdef set get_enclosing_perimeter_verts(self, Vert vert):
        if vert in self.tree_verts:
            raise ValueError("Starting vert is already in tree.")
        if vert.role != VertRole.DUAL:
            raise ValueError("Starting vert role does not match graph_selection.")
            
        cdef set perim = set()
        cdef set visited = set()
        cdef object queue = deque([vert])
        cdef Vert current
        cdef QuadEdge edge
        cdef Vert neighbor
        
        while len(queue) > 0:
            current = queue.popleft()
            visited.add(current)
            
            for edge in current.cc_edges:
                if edge in self.tree_edges:
                    raise RuntimeError("Flood fill leaked into tree edge.")
                neighbor, _ = edge.get_dual_dest_from(current)
                if neighbor not in visited and neighbor not in queue and neighbor not in self.tree_verts:
                    queue.append(neighbor)
                if neighbor in self.tree_verts:
                    perim.add(neighbor)
        return perim

    cpdef set get_enclosed_primal_verts(self, list perim):
        cdef set perim_set = set()
        cdef QuadEdge pe_edge
        for pe_edge, _ in perim:
            perim_set.add(pe_edge)
            
        cdef Vert start_primal = None
        cdef QuadEdge next_edge = (<QuadEdge>perim[0][0]).dual_AB_cc_next
        cdef tuple src_primals = (<QuadEdge>perim[0][0]).get_primal_vert_pair(EdgeDir.AB)
        cdef tuple next_primals = next_edge.get_primal_vert_pair(EdgeDir.AB)
        
        if src_primals[0] in next_primals:
            start_primal = src_primals[0]
        elif src_primals[1] in next_primals:
            start_primal = src_primals[1]
            
        if start_primal is None:
            raise ValueError("Could not find starting primal vert inside perimeter.")
            
        cdef set enclosed = set()
        cdef set visited = set()
        cdef object queue = deque([start_primal])
        cdef Vert current
        cdef QuadEdge edge
        cdef Vert dest
        
        while len(queue) > 0:
            current = queue.popleft()
            visited.add(current)
            enclosed.add(current)
            
            for edge in current.cc_edges:
                if edge not in perim_set:
                    dest, _ = edge.get_primal_dest_from(current)
                    if dest not in visited and dest not in queue:
                        queue.append(dest)
        return enclosed

    cpdef tuple loop_erased_random_walk_from(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region):
        if self.unwind_loops:
            return self._loop_erased_random_walk_from_unwind(vert, existing_walk_edges, existing_consumed_verts, region)
        else:
            return self._loop_erased_random_walk_from_lazy(vert, existing_walk_edges, existing_consumed_verts, region)

    cpdef tuple _loop_erased_random_walk_from_lazy(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region):
        assert vert not in self.tree_verts
        if vert.role >= 0: # Not DUAL or DUAL_EXTERIOR
             raise ValueError("Starting vert role does not match graph_selection.")
             
        cdef QuadEdge current_edge
        cdef Vert current_vert = vert
        cdef dict next_step = {}
        cdef Vert termination_vert = None
        cdef bint escaped_existing_walk = False
        
        cdef Vert next_vert
        cdef Vert end_a, end_b
        
        # Animation variables
        cdef list temp_walk_edges
        cdef Vert temp_curr
        cdef bint temp_escaped
        cdef int trace_steps
        cdef int max_steps
        cdef QuadEdge edge
        cdef Vert next_v
        cdef set perimeter_edges
        cdef list walk_tuples, tree_tuples
        
        while True:
            # Pick next edge
            current_edge = current_vert.cc_edges[random.randint(0, len(current_vert.cc_edges)-1)]
            next_vert, _ = current_edge.get_dual_dest_from(current_vert)
            
            # Record step
            next_step[current_vert] = current_edge
            
            # Check termination
            if next_vert in self.tree_verts:
                termination_vert = next_vert
                break
            
            if escaped_existing_walk and next_vert in existing_consumed_verts:
                termination_vert = next_vert
                break
                
            if next_vert.role == VertRole.DUAL_EXTERIOR:
                termination_vert = next_vert
                break
                
            if ((current_edge.dual_AB is not None and current_edge.dual_AB.role == VertRole.DUAL_EXTERIOR) or
                (current_edge.dual_BA is not None and current_edge.dual_BA.role == VertRole.DUAL_EXTERIOR)):
                termination_vert = next_vert
                break
                
            # Update escaped state
            if not escaped_existing_walk:
                if current_vert not in existing_consumed_verts and next_vert not in existing_consumed_verts:
                    escaped_existing_walk = True
            
            # Move to next vert
            current_vert = next_vert
            
            # Animation tracking
            if self.animating:
                temp_walk_edges = []
                temp_curr = vert
                temp_escaped = False
                trace_steps = 0
                max_steps = len(next_step) + 2
                
                while temp_curr != current_vert and trace_steps < max_steps:
                    if temp_curr not in next_step:
                        break
                    edge = next_step[temp_curr]
                    next_v, _ = edge.get_dual_dest_from(temp_curr)
                    
                    if not temp_escaped:
                        if temp_curr not in existing_consumed_verts and next_v not in existing_consumed_verts:
                            temp_escaped = True
                    
                    if temp_escaped:
                        temp_walk_edges.append(edge)
                        
                    temp_curr = next_v
                    trace_steps += 1
                    
                perimeter_edges = set()
                perimeter_edges.update(region.dual_perimeter_edges)
                walk_tuples = [
                    (edge, VertRole.DUAL, EdgeDir.AB, 0 if edge not in perimeter_edges else 2)
                    for edge in temp_walk_edges + existing_walk_edges
                ]
                tree_tuples = [
                    (edge, VertRole.DUAL, EdgeDir.AB, 0 if edge not in perimeter_edges else 2)
                    for edge in self.tree_edges
                ]
                self.animation_tracks[0].append(walk_tuples + tree_tuples)
                self.animation_tracks[0][-1].extend([
                    (edge, VertRole.DUAL, EdgeDir.AB, 1)
                    for edge in self.bridge_edges
                ])

        # Reconstruct path
        cdef list walk_edges = []
        cdef set consumed_verts = set()
        cdef Vert curr = vert
        cdef bint reconstruct_escaped = False
        
        while curr != termination_vert:
            edge = next_step[curr]
            next_v, _ = edge.get_dual_dest_from(curr)
            
            if not reconstruct_escaped:
                if curr not in existing_consumed_verts and next_v not in existing_consumed_verts:
                    reconstruct_escaped = True
            
            if reconstruct_escaped:
                walk_edges.append(edge)
                consumed_verts.add(curr)
                consumed_verts.add(next_v)
            
            curr = next_v
            
        return walk_edges, consumed_verts, termination_vert

    cpdef tuple _loop_erased_random_walk_from_unwind(self, Vert vert, list existing_walk_edges, set existing_consumed_verts, Region region):
        assert vert not in self.tree_verts
        if vert.role >= 0: # Not DUAL or DUAL_EXTERIOR
             raise ValueError("Starting vert role does not match graph_selection.")
             
        cdef QuadEdge current_edge
        cdef Vert current_vert
        cdef set consumed_verts = set()
        cdef Vert termination_vert = None
        cdef list walk_edges = []
        cdef bint escaped_existing_walk = False
        
        current_edge = vert.cc_edges[random.randint(0, len(vert.cc_edges)-1)]
        current_vert, _ = current_edge.get_dual_dest_from(vert)
        
        cdef Vert rewind_vert
        cdef QuadEdge last_edge
        cdef Vert end_a, end_b
        cdef set perimeter_edges
        cdef list walk_tuples, tree_tuples
        
        while True:
            if current_vert in consumed_verts:
                rewind_vert, _ = current_edge.get_dual_dest_from(current_vert)
                while rewind_vert != current_vert:
                    consumed_verts.remove(rewind_vert)
                    last_edge = walk_edges.pop()
                    rewind_vert, _ = last_edge.get_dual_dest_from(rewind_vert)
            else:
                if escaped_existing_walk:
                    consumed_verts.add(current_vert)
                    walk_edges.append(current_edge)
                
                end_a, end_b = current_edge.get_dual_vert_pair(EdgeDir.AB)
                if not escaped_existing_walk and end_a not in existing_consumed_verts and end_b not in existing_consumed_verts:
                    escaped_existing_walk = True
                    consumed_verts.add(end_a)
                    consumed_verts.add(end_b)
                    walk_edges.append(current_edge)
                    
            if self.animating:
                perimeter_edges = set()
                perimeter_edges.update(region.dual_perimeter_edges)
                walk_tuples = [
                    (edge, VertRole.DUAL, EdgeDir.AB, 0 if edge not in perimeter_edges else 2)
                    for edge in walk_edges + existing_walk_edges
                ]
                tree_tuples = [
                    (edge, VertRole.DUAL, EdgeDir.AB, 0 if edge not in perimeter_edges else 2)
                    for edge in self.tree_edges
                ]
                self.animation_tracks[0].append(walk_tuples + tree_tuples)
                self.animation_tracks[0][-1].extend([
                    (edge, VertRole.DUAL, EdgeDir.AB, 1)
                    for edge in self.bridge_edges
                ])
                
            if (current_vert in self.tree_verts or 
                (current_vert in existing_consumed_verts and escaped_existing_walk) or
                current_vert.role == VertRole.DUAL_EXTERIOR or
                ((current_edge.dual_AB is not None and current_edge.dual_AB.role == VertRole.DUAL_EXTERIOR) or
                 (current_edge.dual_BA is not None and current_edge.dual_BA.role == VertRole.DUAL_EXTERIOR))):
                termination_vert = current_vert
                break
                
            current_edge = current_vert.cc_edges[random.randint(0, len(current_vert.cc_edges)-1)]
            current_vert, _ = current_edge.get_dual_dest_from(current_vert)
            
        return walk_edges, consumed_verts, termination_vert

    cpdef QuadEdge commit_walk(self, list walk_edges, set consumed_verts):
        cdef QuadEdge edge
        cdef Vert vert
        cdef Vert dest
        
        for edge in walk_edges:
            self.tree_edges.add(edge)
        for vert in consumed_verts:
            self.tree_verts.add(vert)
            
        for vert in consumed_verts:
            for edge in vert.cc_edges:
                if edge in self.tree_edges:
                    continue
                dest, _ = edge.get_dual_dest_from(vert)
                if dest in self.tree_verts:
                    self.bridge_edges.add(edge)
                    
        if self.animating:
            self.animation_tracks[0][-1].extend([
                (edge, VertRole.DUAL, EdgeDir.AB, 1)
                for edge in self.bridge_edges
            ])
            
        return walk_edges[0]

    cpdef list develop_region(self, Region src_region, QuadEdge start_edge, EdgeDir start_dir):
        cdef list loop = self.traverse_clockwise_loop(start_edge, start_dir, self.bridge_edges)
        cdef int region_weight = self.graph.count_primal_verts_within_perim(loop)
        cdef Region root_region = Region(region_weight, loop)
        self.region_tree.add_region(root_region)
        
        cdef list new_regions = [root_region]
        cdef QuadEdge edge
        cdef EdgeDir dir
        cdef list child_regions
        cdef Region neighbor_region
        cdef Edge new_region_edge
        cdef Edge old_region_edge
        cdef Region other_region
        
        for edge, dir in loop:
            if edge in self.bridge_edges and edge not in src_region.bridge_set and edge != start_edge:
                child_regions = self.develop_region(src_region, edge, reverse_dir(dir))
                new_regions.extend(child_regions)
                neighbor_region = child_regions[0]
                
                new_region_edge = Edge()
                new_region_edge.end_A = root_region
                new_region_edge.end_B = neighbor_region
                new_region_edge.twin_graph_edge = edge
                
                self.region_tree.add_edge(new_region_edge)
                new_region_edge.calculate_weight_differential_from_dir(EdgeDir.AB)
                
            if edge in self.bridge_edges and edge in src_region.bridge_set:
                old_region_edge = src_region.bridge_to_region_edge_map[edge]
                other_region, _ = old_region_edge.get_dest_from(src_region)
                self.region_tree.remove_edge(old_region_edge)
                
                new_region_edge = Edge()
                new_region_edge.end_A = root_region
                new_region_edge.end_B = other_region
                new_region_edge.twin_graph_edge = edge
                
                self.region_tree.add_edge(new_region_edge)
                new_region_edge.calculate_weight_differential_from_dir(EdgeDir.AB)
                
        return new_regions

    cpdef list traverse_clockwise_loop(self, QuadEdge start_edge, EdgeDir edge_dir, set bridge_boundaries):
        cdef list visited_edges = []
        cdef Vert current_vert
        cdef QuadEdge current_edge = start_edge
        cdef EdgeDir current_dir = edge_dir
        
        current_vert, _ = start_edge.get_dual_vert_pair(edge_dir)
        cdef Vert start_vert = current_vert
        
        cdef QuadEdge next_edge
        cdef EdgeDir next_dir
        
        while True:
            current_vert, current_dir = current_edge.get_dual_dest_from(current_vert)
            visited_edges.append((current_edge, current_dir))
            
            next_edge = current_edge
            next_dir = current_dir
            while True:
                if current_vert == next_edge.dual_BA:
                    next_edge = next_edge.dual_BA_cc_next
                else:
                    next_edge = next_edge.dual_AB_cc_next
                _, next_dir = next_edge.get_dual_dest_from(current_vert)
                
                if next_edge in self.tree_edges or next_edge in bridge_boundaries:
                    break
            current_edge = next_edge
            current_dir = next_dir
            
            if current_edge == start_edge and current_vert == start_vert:
                break
                
        return visited_edges
