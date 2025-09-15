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
        red = (255,100,100)
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

        # Loop
        exitCondition = False
        while not exitCondition:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    exitCondition = True

            # Screen Reset
            screen.fill(GraphVis.DrawColors.black)

            # Mouse
            mx, my = pg.mouse.get_pos()
            w_mp = Point(mx, my)
            g_mp = self.point_window_to_graph(w_mp)
            # pg.draw.circle(screen, GraphVis.DrawColors.blue, w_mp.tuple(), 5)

            # Draw edges
            for vert in self.graph.primalVerts:
                for edge in vert.cc_edges:
                    dest, dir = edge.getPrimalDestFrom(vert)
                    if dir == TwinGraph.EdgeDir.AB: # Only draw one side of the graph to avoid redundancy
                        src = self.point_graph_to_window(vert.point)
                        dest = self.point_graph_to_window(dest.point)
                        pg.draw.aaline(screen, GraphVis.DrawColors.white, src.tuple(), dest.tuple(), blend=10)

            # Draw Close Vertices
            closest_vert = graph.get_closest_vert(g_mp, TwinGraph.VertRole.PRIMAL)
            pg.draw.circle(screen, GraphVis.DrawColors.white, self.point_graph_to_window(closest_vert.point).tuple(), 5)

            # Push to screen
            pg.display.flip()

        # Quit Pygame
        pg.quit()

    def point_graph_to_window(self, point: Point) -> Point:
        transformed_coords = self.graph_to_window(point.x, point.y) 
        return Point(transformed_coords[0], transformed_coords[1])
    def point_window_to_graph(self, point: Point) -> Point:
        transformed_coords = self.window_to_graph(point.x, point.y) 
        return Point(transformed_coords[0], transformed_coords[1])

    def graph_to_window(self, graph_x: float, graph_y: float) -> (float, float):
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
    
    def window_to_graph(self, window_x: float, window_y: float) -> (float, float):
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

    