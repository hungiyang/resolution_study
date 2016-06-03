__all__ = ['count_pairs_JKindex', 'point2JKindex']

import numpy as np
from fast3tree import fast3tree

def point2JKindex(loc_array, box_size = 250., n = 50):
    dL = box_size/float(n)
    index = np.floor_divide(loc_array, dL).astype(int)
    return np.ravel_multi_index(index.T, (n,n,n))

def count_pairs_JKindex(points, rbins, box_size):
    pairs = [[] for r in rbins]
    with fast3tree(points) as tree:
        tree.set_boundaries(0, box_size)
        for pidx, p in enumerate(points):
            pJKidx = point2JKindex(p, box_size = box_size)
            for i, r in enumerate(rbins):
                idx, pos = tree.query_radius(p,r, periodic=True, output='both')
                pos = pos[idx < pidx]
                pairs[i].extend([pJKidx, j] for j in point2JKindex(pos, box_size = box_size))
    pairs = np.array([np.array(i) for i in pairs])
    return np.array(pairs)


#def count_pairs_JKindex_old(points, rbins, box_size):
#    pairs = np.zeros(len(rbins), dtype=int)
#    with fast3tree(points) as tree:
#        tree.set_boundaries(0, box_size)
#        pairs = []
#        for i, r in enumerate(rbins):
#            loc = []
#            for p in points:
#                pindex = point2JKindex(p)
#                #save the distinct pairs
#                ppairs = np.array([np.sort([pindex, point2JKindex(pfound)]) for pfound in tree.query_radius(p,r, periodic=True, output='pos') if np.sum(np.abs(pfound - p)) > 0.0]) 
#                if len(ppairs) :
#                    loc.append(ppairs)
#            #Since every pair is counted twice
#            if len(loc) > 0:
#                loc = np.sort(np.concatenate(loc), axis = 0)[::2]
#            else:
#                loc = np.array([])
#            pairs.append(loc)
#    return np.array(pairs)