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
        displaying_labels = False
        displaying_region_tree = False
        displaying_perimeters = False
        selected_vert: Optional[TwinGraph.Vert] = None

        # Animation
        animation_frame = 0
        animation_paused = False
        animation_deck_idx = 0
        animation_deck = self.graph.animation_tracks + self.graph_nav_dual.animation_tracks
        for track in animation_deck:
            assert len(track) > 0, "Animation deck contains empty track."

        # Create label surfaces
        font = pg.font.SysFont(None, 13)
        # Create a transparent surface for all labels
        primal_labels_surface = pg.Surface(GraphVis.window_size, pg.SRCALPHA)
        dual_labels_surface = pg.Surface(GraphVis.window_size, pg.SRCALPHA)
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
                animation_frame = animation_frame-5
            if keys[pg.K_RIGHT]:
                animation_frame = animation_frame+5
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
                    src, dest = self.get_exterior_adjusted_point(edge.primal_A, edge.primal_B)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                    pg.draw.aaline(screen, GraphVis.DrawColors.white, src.tuple(), dest.tuple(), blend=10)
                if displaying_dual and edge.dual_AB is not None and edge.dual_BA is not None: # Check dual exists and should be drawn
                    if displaying_dual_external or ( # Check whether dual edge is external and should be drawn
                        edge.dual_AB.role != TwinGraph.VertRole.DUAL_EXTERIOR and
                        edge.dual_BA.role != TwinGraph.VertRole.DUAL_EXTERIOR
                    ):
                            src, dest = self.get_exterior_adjusted_point(edge.dual_AB, edge.dual_BA)
                            src, dest = self.point_pair_graph_to_window(src, dest)
                            pg.draw.aaline(screen, GraphVis.DrawColors.blue, src.tuple(), dest.tuple(), blend=10)

            # Draw Animation
            for edge, role, _, tag in animation_deck[animation_deck_idx][animation_frame]:
                endA, endB = edge.get_dual_vert_pair(TwinGraph.EdgeDir.AB)
                if endA.role == TwinGraph.VertRole.DUAL_EXTERIOR or endB.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                    if not displaying_dual_external:
                        continue

                if role == TwinGraph.VertRole.PRIMAL:
                    src, dest = self.get_exterior_adjusted_point(edge.primal_A, edge.primal_B)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                else:
                    if edge.dual_AB is None or edge.dual_BA is None:
                        continue
                    src, dest = self.get_exterior_adjusted_point(edge.dual_AB, edge.dual_BA)
                    src, dest = self.point_pair_graph_to_window(src, dest)
                
                color: Tuple[int,int,int] = GraphVis.DrawColors.white # Placeholder
                if tag == 0:
                    color = GraphVis.DrawColors.red
                elif tag == 1:
                    color = GraphVis.DrawColors.green
                elif tag == 2:
                    color = GraphVis.DrawColors.purple
                pg.draw.line(screen, color, src.tuple(), dest.tuple(), 5)

            # Draw vertex labels
            if displaying_labels:
                if displaying_primal and not displaying_dual:
                    screen.blit(primal_labels_surface, (0,0))
                if displaying_dual and not displaying_primal:
                    screen.blit(dual_labels_surface, (0,0))

            # Draw Region Tree
            if displaying_region_tree:
                for edge in self.region_tree.edges:
                    src, dest = self.point_pair_graph_to_window(edge.end_A.point, edge.end_B.point)
                    pg.draw.aaline(screen, GraphVis.DrawColors.green, src.tuple(), dest.tuple(), blend=10)
                for region in self.region_tree.regions:
                    region_pos = self.point_graph_to_window(region.point)
                    label_surface = font.render(str(region.id_str), True, GraphVis.DrawColors.white)
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
                            src, dest = self.get_exterior_adjusted_point(*edge.get_dual_vert_pair(dir))
                            src, dest = self.point_pair_graph_to_window(src, dest)
                            pg.draw.line(screen, GraphVis.DrawColors.pink, src.tuple(), dest.tuple(), 5)

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
    
    def get_exterior_adjusted_point(self, vertA: TwinGraph.Vert, vertB: TwinGraph.Vert) -> Tuple[Point, Point]:
        if vertA.role == TwinGraph.VertRole.DUAL_EXTERIOR or vertB.role == TwinGraph.VertRole.DUAL_EXTERIOR:
            # Identify the exterior and non-exterior vertices
            if vertA.role == TwinGraph.VertRole.DUAL_EXTERIOR:
                exterior_vert = vertA
                non_exterior_vert = vertB
            else:
                exterior_vert = vertB
                non_exterior_vert = vertA

            # Center of the graph region
            center_x = (self.graph.upperXBound + self.graph.lowerXBound) / 2
            center_y = (self.graph.upperYBound + self.graph.lowerYBound) / 2
            center = Point(center_x, center_y)

            # Vector from center to non-exterior point
            dx = non_exterior_vert.point.x - center.x
            dy = non_exterior_vert.point.y - center.y

            # Angle from center to non-exterior point
            angle = math.atan2(dy, dx)

            # Extend out at that angle beyond the window
            # Use a distance larger than the graph bounds
            extend_dist = max(
                self.graph.upperXBound - self.graph.lowerXBound,
                self.graph.upperYBound - self.graph.lowerYBound
            ) * 2

            dest_x = non_exterior_vert.point.x + extend_dist * math.cos(angle)
            dest_y = non_exterior_vert.point.y + extend_dist * math.sin(angle)
            dest = Point(dest_x, dest_y)

            return non_exterior_vert.point, dest
        else:
            return vertA.point, vertB.point

    