import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

import numpy as np

from Euclid import Point
from typing import List, Tuple

def generate_grid_graph(width: int, height: int) -> Tuple[List[Point], List[int], List[Tuple[int, int]]]:
    w, h = int(width), int(height)
    n = w * h

    # Points
    xs = np.tile(np.arange(w, dtype=int), h)
    ys = np.repeat(np.arange(h, dtype=int), w)
    points = [Point(int(x), int(y)) for x, y in zip(xs, ys)]

    weights = [1] * n

    # Horizontal edges
    if w > 1 and h > 0:
        start_h = (np.arange(h, dtype=int)[:, None] * w + np.arange(w - 1, dtype=int)[None, :]).ravel()
        end_h = start_h + 1
        horiz = np.stack((start_h, end_h), axis=1)
    else:
        horiz = np.empty((0, 2), dtype=int)

    # Vertical edges
    if h > 1 and w > 0:
        start_v = (np.arange(h - 1, dtype=int)[:, None] * w + np.arange(w, dtype=int)[None, :]).ravel()
        end_v = start_v + w
        vert = np.stack((start_v, end_v), axis=1)
    else:
        vert = np.empty((0, 2), dtype=int)

    edges = np.vstack((horiz, vert))
    edgeIdxs = [(int(a), int(b)) for a, b in edges]

    return points, weights, edgeIdxs