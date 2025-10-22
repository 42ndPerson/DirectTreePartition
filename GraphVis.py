from typing import List, Tuple, Set, Optional

import pygame as pg

from TwinGraph import *
from Euclid import *
from GraphNav import *
from RegionTree import *
import math

class GraphVis:
    content_padding = 50
    window_size = (800,800)

    # Colors
    class DrawColors:
        black = (0,0,0)
        white = (255,255,255)
        red = (255,50,100)
        green = (100,255,100)
        blue = (100,100,255)
        purple = (200,100,255)
        pink = (255,50,255)
        orange = (255,70,0)

    def __init__(self, graph: TwinGraph, graph_nav_dual: GraphNav, region_tree: RegionTree) -> None:
        # Graph
        self.graph = graph

        # GraphNav
        self.graph_nav_dual = graph_nav_dual

        # RegionTree
        self.region_tree = region_tree

        # Init
        pg.init()

        # Screen
        screen = pg.display.set_mode(GraphVis.window_size)
        pg.display.set_caption("Graph Viewer")

        # Loop State
        loop_idx = 0
        displaying_primal = True
        displaying_dual = False
        displaying_dual_external = True
        displaying_dual_annotations = False
        displaying_labels = False
        displaying_region_tree = False
        displaying_perimeters = False
        displaying_central_region = False
        selected_vert: Optional[TwinGraph.Vert] = None

        # Animation
        animation_frame = 0
        animation_paused = False
        animation_deck_idx = 4
        animation_deck = self.graph.animation_tracks + self.graph_nav_dual.animation_tracks
        for track in animation_deck:
            assert len(track) > 0, "Animation deck contains empty track."

        # Create label surfaces
        font = pg.font.SysFont(None, 13)
        # Create a transparent surface for all labels
        primal_labels_surface = pg.Surface(GraphVis.window_size, pg.SRCALPHA)
        dual_labels_surface = pg.Surface(GraphVis.window_size, pg.SRCALPHA)
        dual_annotations_labels_surface = pg.Surface(GraphVis.window_size, pg.SRCALPHA)
        for vert in self.graph.primalVerts:
            vert_pos = self.point_graph_to_window(vert.point)
            label_pos = (vert_pos.x + 8, vert_pos.y - 12)
            label_surface = font.render(str(vert.id_str), True, GraphVis.DrawColors.purple)
            primal_labels_surface.blit(label_surface, label_pos)
        for vert in self.graph.dualVerts:
            vert_pos = self.point_graph_to_window(vert.point)
            label_pos = (vert_pos.x + 8, vert_pos.y - 12)
            label_surface = font.render(str(vert.id_str), True, GraphVis.DrawColors.green)
            dual_labels_surface.blit(label_surface, label_pos)
        for edge in self.graph.edges:
            if edge.dual_AB is not None and edge.dual_BA is not None: # Check dual exists and should be drawn
                src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.DUAL)
                src, dest = self.point_pair_graph_to_window(src, dest)
        
                # Put text label for dual edge annotation
                # Compute midpoint and scale toward origin if outside screen bounds
                buffer = 20
                mid_x = (src.x + dest.x) / 2
                mid_y = (src.y + dest.y) / 2
                max_x = GraphVis.window_size[0] - buffer
                max_y = GraphVis.window_size[1] - buffer
                min_x = buffer
                min_y = buffer
                origin_x = GraphVis.window_size[0] / 2
                origin_y = GraphVis.window_size[1] / 2

                # If midpoint is outside bounds, scale toward origin
                def scale_toward_origin(x, y, min_x, max_x, min_y, max_y, origin_x, origin_y):
                    if min_x <= x <= max_x and min_y <= y <= max_y:
                        return x, y
                    dx = x - origin_x
                    dy = y - origin_y
                    scale_x = 1.0
                    scale_y = 1.0
                    if x < min_x:
                        scale_x = (min_x - origin_x) / dx if dx != 0 else 1.0
                    elif x > max_x:
                        scale_x = (max_x - origin_x) / dx if dx != 0 else 1.0
                    if y < min_y:
                        scale_y = (min_y - origin_y) / dy if dy != 0 else 1.0
                    elif y > max_y:
                        scale_y = (max_y - origin_y) / dy if dy != 0 else 1.0
                    scale = min(scale_x, scale_y)
                    return origin_x + dx * scale, origin_y + dy * scale

                mid_x, mid_y = scale_toward_origin(mid_x, mid_y, min_x, max_x, min_y, max_y, origin_x, origin_y)
                mid_point = Point(mid_x, mid_y)
                if edge.dual_AB_annotation is not None:
                    label_surface = font.render(str(edge.dual_AB_annotation), True, GraphVis.DrawColors.orange)
                    dual_annotations_labels_surface.blit(label_surface, (mid_point.x, mid_point.y)) 
                else:
                    label_surface = font.render("NA", True, GraphVis.DrawColors.orange)
                    dual_annotations_labels_surface.blit(label_surface, (mid_point.x, mid_point.y)) 

        # Loop
        exitCondition = False
        while not exitCondition:
            # Key / System Input
            # Quitting / Toggles
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    exitCondition = True
                if event.type == pg.KEYDOWN:
                    if event.key == pg.K_SPACE:
                        animation_paused = not animation_paused
                    if event.key == pg.K_p:
                        displaying_primal = not displaying_primal
                    if event.key == pg.K_d:
                        displaying_dual = not displaying_dual
                    if event.key == pg.K_e:
                        displaying_dual_external = not displaying_dual_external
                    if event.key == pg.K_l:
                        displaying_labels = not displaying_labels
                    if event.key == pg.K_r:
                        displaying_region_tree = not displaying_region_tree
                    if event.key == pg.K_b:
                        displaying_perimeters = not displaying_perimeters
                    if event.key == pg.K_a:
                        displaying_dual_annotations = not displaying_dual_annotations
                    if event.key == pg.K_c:
                        displaying_central_region = not displaying_central_region
                    if event.key == pg.K_q:
                        self.graph_nav_dual.run_two_split_attempt()
                    if event.key == pg.K_m:
                        # Reset system
                        self.region_tree = RegionTree(self.graph)
                        self.graph_nav_dual = GraphNav(self.graph, self.region_tree)
                        animation_deck = self.graph.animation_tracks + self.graph_nav_dual.animation_tracks
                    if pg.K_1 <= event.key <= pg.K_9: # Select animation track
                        idx = event.key - pg.K_1
                        if idx < len(animation_deck):
                            animation_deck_idx = idx
                            animation_frame = 0
                if event.type == pg.MOUSEBUTTONDOWN:
                    if not displaying_perimeters:
                        # Handle mouse or trackpad clicks
                        selected_vert = self.graph.get_closest_vert(
                            self.point_window_to_graph(Point(*pg.mouse.get_pos())),
                            TwinGraph.VertRole.DUAL if displaying_dual else TwinGraph.VertRole.PRIMAL
                        )
                        if selected_vert is not None and selected_vert.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                            selected_vert = None # Can't start from exterior dual verts
                        if selected_vert is not None and selected_vert.role == TwinGraph.VertRole.DUAL:
                            # self.graph_nav_dual.loop_erased_random_walk_from(selected_vert)
                            selected_vert = self.graph.get_closest_vert(
                                self.point_window_to_graph(Point(*pg.mouse.get_pos())),
                                TwinGraph.VertRole.DUAL if displaying_dual else TwinGraph.VertRole.PRIMAL
                            )
                            if selected_vert is not None and selected_vert.role == TwinGraph.VertRole.DUAL:
                                self.graph_nav_dual.walk_division_from(None, selected_vert)

            # Sustained Actions
            keys = pg.key.get_pressed()
            if keys[pg.K_LEFT]:
                animation_frame = animation_frame-3
            if keys[pg.K_RIGHT]:
                animation_frame = animation_frame+3
            if animation_frame < 0:
                animation_frame = len(animation_deck[animation_deck_idx]) + animation_frame
            if animation_frame >= len(animation_deck[animation_deck_idx]):
                animation_frame = 0

            # Screen Reset
            screen.fill(GraphVis.DrawColors.black)

            # Mouse
            mx, my = pg.mouse.get_pos()
            w_mp = Point(mx, my)
            g_mp = self.point_window_to_graph(w_mp)

            # Draw edges
            for edge in self.graph.edges:
                if displaying_primal and edge in self.graph.edges.difference(set(animation_deck[animation_deck_idx][animation_frame])):
                    src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.PRIMAL)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                    pg.draw.aaline(screen, GraphVis.DrawColors.white, src.tuple(), dest.tuple(), blend=10)
                if displaying_dual and edge.dual_AB is not None and edge.dual_BA is not None: # Check dual exists and should be drawn
                    if displaying_dual_external or ( # Check whether dual edge is external and should be drawn
                        edge.dual_AB.role != TwinGraph.VertRole.DUAL_EXTERIOR and
                        edge.dual_BA.role != TwinGraph.VertRole.DUAL_EXTERIOR
                    ):
                            src, dest, dir = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.DUAL)
                            src, dest = self.point_pair_graph_to_window(src, dest)
                            pg.draw.aaline(screen, GraphVis.DrawColors.blue, src.tuple(), dest.tuple(), blend=10)

                            if displaying_dual_annotations:
                                # Draw arrow on dual edge indicating AB direction
                                # If one vertex is exterior, draw arrow at the non-exterior vertex pointing outward
                                arrow_length = 6
                                arrow_width = 3
                                if edge.dual_AB.role == TwinGraph.VertRole.DUAL_EXTERIOR or edge.dual_BA.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                                    # Find non-exterior vertex and direction
                                    adjusted_src = src
                                    adjusted_dest = dest
                                    # If src is outside the screen area, swap adjusted_src and adjusted_dest
                                    buffer = 20
                                    if not (buffer <= adjusted_src.x <= GraphVis.window_size[0] - buffer and buffer <= adjusted_src.y <= GraphVis.window_size[1] - buffer):
                                        adjusted_src, adjusted_dest = adjusted_dest, adjusted_src
                                        dir = TwinGraph.EdgeDir.BA if dir == TwinGraph.EdgeDir.AB else TwinGraph.EdgeDir.AB

                                    dx = adjusted_src.x - adjusted_dest.x
                                    dy = adjusted_src.y - adjusted_dest.y
                                    dx /= 20
                                    dy /= 20
                                    length = math.hypot(dx, dy)
                                    if length > 0:
                                        tip_x = adjusted_src.x - dx * 0.18
                                        tip_y = adjusted_src.y - dy * 0.18
                                        if dir == TwinGraph.EdgeDir.AB:
                                            base_x = tip_x + dx / length * arrow_length
                                            base_y = tip_y + dy / length * arrow_length
                                        else:
                                            base_x = tip_x - dx / length * arrow_length
                                            base_y = tip_y - dy / length * arrow_length
                                        perp_x = -dy / length * arrow_width
                                        perp_y = dx / length * arrow_width
                                        left_x = base_x + perp_x
                                        left_y = base_y + perp_y
                                        right_x = base_x - perp_x
                                        right_y = base_y - perp_y
                                        pg.draw.aaline(screen, GraphVis.DrawColors.orange, (tip_x, tip_y), (left_x, left_y), blend=10)
                                        pg.draw.aaline(screen, GraphVis.DrawColors.orange, (tip_x, tip_y), (right_x, right_y), blend=10)
                                else:
                                    # Standard arrow in AB direction
                                    dx = dest.x - src.x
                                    dy = dest.y - src.y
                                    length = math.hypot(dx, dy)
                                    if length > 0:
                                        tip_x = dest.x - dx * 0.18
                                        tip_y = dest.y - dy * 0.18
                                        base_x = tip_x - dx / length * arrow_length
                                        base_y = tip_y - dy / length * arrow_length
                                        perp_x = -dy / length * arrow_width
                                        perp_y = dx / length * arrow_width
                                        left_x = base_x + perp_x
                                        left_y = base_y + perp_y
                                        right_x = base_x - perp_x
                                        right_y = base_y - perp_y
                                        pg.draw.aaline(screen, GraphVis.DrawColors.orange, (tip_x, tip_y), (left_x, left_y), blend=10)
                                        pg.draw.aaline(screen, GraphVis.DrawColors.orange, (tip_x, tip_y), (right_x, right_y), blend=10)

            # Draw Animation
            for edge, role, _, tag in animation_deck[animation_deck_idx][animation_frame]:
                endA, endB = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
                if endA.role == TwinGraph.VertRole.DUAL_EXTERIOR or endB.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                    if not displaying_dual_external:
                        continue

                if role == TwinGraph.VertRole.PRIMAL:
                    src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.PRIMAL)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                else:
                    src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.DUAL)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                
                color: Tuple[int,int,int] = GraphVis.DrawColors.white # Placeholder
                if tag == 0:
                    color = GraphVis.DrawColors.red
                elif tag == 1:
                    color = GraphVis.DrawColors.green
                elif tag == 2:
                    color = GraphVis.DrawColors.purple
                pg.draw.line(screen, color, src.tuple(), dest.tuple(), 5)

            # Highlight central region
            if displaying_central_region and self.region_tree.central_region is not None:
                for edge, dir in self.region_tree.central_region.dual_perimeter:
                    src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.DUAL)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                    pg.draw.line(screen, GraphVis.DrawColors.pink, src.tuple(), dest.tuple(), 5)

            # Draw vertex labels
            if displaying_labels:
                if displaying_primal and not displaying_dual:
                    screen.blit(primal_labels_surface, (0,0))
                if displaying_dual and not displaying_primal:
                    screen.blit(dual_labels_surface, (0,0))
            # Draw dual edge annotation labels
            if displaying_dual and displaying_dual_annotations:
                screen.blit(dual_annotations_labels_surface, (0,0))

            # Draw Region Tree
            if displaying_region_tree:
                for edge in self.region_tree.edges:
                    alert = False
                    if edge.twin_graph_edge in map(lambda e: e.twin_graph_edge, self.region_tree.edges - {edge}):
                        alert = True

                    src, dest = self.point_pair_graph_to_window(edge.end_A.point, edge.end_B.point)
                    pg.draw.aaline(screen, GraphVis.DrawColors.green if not alert else GraphVis.DrawColors.orange, src.tuple(), dest.tuple(), blend=10)

                    if edge.ab_weight_differential is not None:
                        mid_x = (src.x + dest.x) / 2
                        mid_y = (src.y + dest.y) / 2
                        label_surface = font.render(f"{edge.ab_weight_differential:.0f}", True, GraphVis.DrawColors.red if not alert else GraphVis.DrawColors.orange)
                        screen.blit(label_surface, (mid_x+2, mid_y+2))
                for region in self.region_tree.regions:
                    region_pos = self.point_graph_to_window(region.point)
                    label_surface = font.render(str(region.id_str) + "---" + str(region.weight), True, GraphVis.DrawColors.white)
                    screen.blit(label_surface, (region_pos.x + 8, region_pos.y - 12))
                # for region in self.region_tree.regions:
                #     pg.draw.circle(screen, GraphVis.DrawColors.red, self.point_graph_to_window(region.point).tuple(), 5)

            # Highlight Selected Vert
            if not displaying_perimeters:
                selected_vert = graph.get_closest_vert(g_mp, TwinGraph.VertRole.DUAL if displaying_dual else TwinGraph.VertRole.PRIMAL)
                if selected_vert is not None:
                    pg.draw.circle(screen, GraphVis.DrawColors.white, self.point_graph_to_window(selected_vert.point).tuple(), 5)
            else:
                selected_region: Optional[RegionTree.Region] = None
                if g_mp is not None:
                    closest_region = None
                    closest_dist = float('inf')
                    for region in self.region_tree.regions:
                        dist = Point.dist(region.point, g_mp)
                        if dist < closest_dist:
                            closest_dist = dist
                            closest_region = region
                    if closest_region is not None and closest_dist < 1.0: # Threshold distance to select region
                        selected_region = closest_region
                if selected_region is not None:
                    if displaying_perimeters:
                        for edge, dir in selected_region.dual_perimeter:
                            src, dest, _ = self.get_exterior_adjusted_point(edge, TwinGraph.VertRole.DUAL)
                            src, dest = self.point_pair_graph_to_window(src, dest)
                            pg.draw.line(screen, GraphVis.DrawColors.pink, src.tuple(), dest.tuple(), 5)
                        for vert in selected_region.get_interior_lining_verts():
                            pg.draw.circle(screen, GraphVis.DrawColors.white, self.point_graph_to_window(vert.point).tuple(), 5)

            # Push to screen
            pg.display.flip()

            # Temporal State Updates
            loop_idx += 1
            if loop_idx%5 == 0 and not animation_paused:
                animation_frame += 1
                if animation_frame >= len(animation_deck[animation_deck_idx]):
                    animation_frame = len(animation_deck[animation_deck_idx]) - 1

        # Quit Pygame
        pg.quit()

    def point_pair_graph_to_window(self, pointA: Point, pointB: Point) -> Tuple[Point, Point]:
        return (self.point_graph_to_window(pointA), self.point_graph_to_window(pointB))
    def point_pair_window_to_graph(self, pointA: Point, pointB: Point) -> Tuple[Point, Point]:
        return (self.point_window_to_graph(pointA), self.point_window_to_graph(pointB))
    def point_graph_to_window(self, point: Point) -> Point:
        transformed_coords = self.graph_to_window(point.x, point.y) 
        return Point(transformed_coords[0], transformed_coords[1])
    def point_window_to_graph(self, point: Point) -> Point:
        transformed_coords = self.window_to_graph(point.x, point.y) 
        return Point(transformed_coords[0], transformed_coords[1])

    def graph_to_window(self, graph_x: float, graph_y: float) -> Tuple[float, float]:
        active_window = (
            GraphVis.window_size[0] - 2*GraphVis.content_padding, 
            GraphVis.window_size[1] - 2*GraphVis.content_padding
        )
        graph_area = (
            self.graph.upperXBound - self.graph.lowerXBound,
            self.graph.upperYBound - self.graph.lowerYBound,
        )

        active_x = active_window[0] * (graph_x - self.graph.lowerXBound) / graph_area[0]
        active_y = active_window[1] * (graph_y - self.graph.lowerYBound) / graph_area[1]

        return (
            GraphVis.content_padding + active_x,
            GraphVis.window_size[1] - GraphVis.content_padding - active_y
        )
    
    def window_to_graph(self, window_x: float, window_y: float) -> Tuple[float, float]:
        active_window = (
            GraphVis.window_size[0] - 2*GraphVis.content_padding, 
            GraphVis.window_size[1] - 2*GraphVis.content_padding
        )
        graph_area = (
            self.graph.upperXBound - self.graph.lowerXBound,
            self.graph.upperYBound - self.graph.lowerYBound,
        )

        active_x = graph_area[0] * (window_x - GraphVis.content_padding) / active_window[0]
        active_y = graph_area[1] * (GraphVis.window_size[1] - GraphVis.content_padding - window_y) / active_window[1]

        return (
            active_x + self.graph.lowerXBound,
            active_y + self.graph.lowerYBound
        )

    def get_exterior_adjusted_point(self, edge: TwinGraph.QuadEdge, role: TwinGraph.VertRole) -> Tuple[Point, Point, TwinGraph.EdgeDir]:
        if role.is_dual():
            vert_A = edge.dual_AB
            vert_B = edge.dual_BA
        else:
            vert_A = edge.primal_A
            vert_B = edge.primal_B

        if vert_A is not None and vert_B is not None:
            if vert_A.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                angle, _ = edge.get_dual_rad_from(vert_B)

                # Project a point out along the angle given by get_primal_rad_from or get_dual_rad_from
                projected_point = Point(
                    vert_B.point.x + 20*math.cos(angle),
                    vert_B.point.y + 20*math.sin(angle)
                )

                return (vert_B.point, projected_point, TwinGraph.EdgeDir.BA)
            
            if vert_B.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                angle, _ = edge.get_dual_rad_from(vert_A)

                # Project a point out along the angle given by get_primal_rad_from or get_dual_rad_from
                projected_point = Point(
                    vert_A.point.x + 20*math.cos(angle),
                    vert_A.point.y + 20*math.sin(angle)
                )

                return (vert_A.point, projected_point, TwinGraph.EdgeDir.AB)
            
            # Neither vertex is exterior, return original points
            return (vert_A.point, vert_B.point, TwinGraph.EdgeDir.AB)
        
        else:
            raise ValueError("Edge does not have both dual vertices defined.")
    