__all__ = ['count_pairs']

import numpy as np
from fast3tree import fast3tree

def count_pairs(points, rbins, box_size):
    pairs = np.zeros(len(rbins), dtype=int)
    with fast3tree(points) as tree:
        tree.set_boundaries(0, box_size)
        for p in points:
            pairs += np.array([tree.query_radius(p, r, periodic=True, \
                    output='c') for r in rbins])
    return (pairs-points.shape[0])/2

