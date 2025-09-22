from typing import List, Tuple, Set

import pygame as pg

from TwinGraph import *
from Euclid import *

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

    def __init__(self, graph: TwinGraph) -> None:
        # Graph
        self.graph = graph

        # Init
        pg.init()

        # Screen
        screen = pg.display.set_mode(GraphVis.window_size)
        pg.display.set_caption("Graph Viewer")

        # Loop State
        loop_idx = 0
        animation_frame = 0
        animation_paused = False
        displaying_primal = True
        displaying_dual = False
        displaying_dual_external = True

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
            # Sustained Actions
            keys = pg.key.get_pressed()
            if keys[pg.K_LEFT]:
                animation_frame = animation_frame-5
            if keys[pg.K_RIGHT]:
                animation_frame = animation_frame+5
            if animation_frame < 0:
                animation_frame = len(self.graph.animation_track) + animation_frame
            if animation_frame >= len(self.graph.animation_track):
                animation_frame = 0

            # Screen Reset
            screen.fill(GraphVis.DrawColors.black)

            # Mouse
            mx, my = pg.mouse.get_pos()
            w_mp = Point(mx, my)
            g_mp = self.point_window_to_graph(w_mp)

            # Draw edges
            for edge in self.graph.edges:
                if displaying_primal and edge in self.graph.edges.difference(set(self.graph.animation_track[animation_frame])):
                    src = self.point_graph_to_window(edge.primal_A.point)
                    dest = self.point_graph_to_window(edge.primal_B.point)
                    pg.draw.aaline(screen, GraphVis.DrawColors.white, src.tuple(), dest.tuple(), blend=10)
                if displaying_dual and edge.dual_AB is not None and edge.dual_BA is not None: # Check dual exists and should be drawn
                    if displaying_dual_external or ( # Check whether dual edge is external and should be drawn
                        edge.dual_AB.role != TwinGraph.VertRole.DUAL_EXTERIOR and
                        edge.dual_BA.role != TwinGraph.VertRole.DUAL_EXTERIOR
                    ):
                            src = self.point_graph_to_window(edge.dual_AB.point)
                            dest = self.point_graph_to_window(edge.dual_BA.point)
                            pg.draw.aaline(screen, GraphVis.DrawColors.blue, src.tuple(), dest.tuple(), blend=10)

            # Test
            for edge, _, _ in self.graph.animation_track[animation_frame]:
                src = self.point_graph_to_window(edge.primal_A.point)
                dest = self.point_graph_to_window(edge.primal_B.point)
                pg.draw.aaline(screen, GraphVis.DrawColors.red, src.tuple(), dest.tuple(), blend=10)

            # Draw Close Vertices
            closest_vert = graph.get_closest_vert(g_mp, TwinGraph.VertRole.PRIMAL)
            if closest_vert is not None:
                pg.draw.circle(screen, GraphVis.DrawColors.white, self.point_graph_to_window(closest_vert.point).tuple(), 5)

            # Push to screen
            pg.display.flip()

            # Temporal State Updates
            loop_idx += 1
            if loop_idx%5 == 0 and not animation_paused:
                animation_frame += 1
                if animation_frame >= len(self.graph.animation_track):
                    animation_frame = 0

        # Quit Pygame
        pg.quit()

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

    