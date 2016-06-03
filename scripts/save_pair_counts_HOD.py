from argparse import ArgumentParser
import numpy as np
from fast3tree import fast3tree

parser = ArgumentParser()
parser.add_argument('case')
args = parser.parse_args()
case = args.case

def count_pairs(points, rbins, box_size):
    pairs = np.zeros(len(rbins), dtype=int)
    with fast3tree(points) as tree:
        tree.set_boundaries(0, box_size)
        for p in points:
            pairs += np.array([tree.query_radius(p, r, periodic=True, \
                    output='c') for r in rbins])
    return (pairs-points.shape[0])/2


#settings
rbins = np.logspace(-1.3, 1.3, 27)
m_log_list = np.linspace(14.5, 10.5, 9)

box_size = case[1:4]
hlist_path ='/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case)

#run
halos = np.load(hlist_path)
halos = halos[halos['upid']==-1]
proxy = 'mvir'
s = halos[proxy].argsort()
k = np.searchsorted(halos[proxy], 10.**m_log_list, sorter=s)

halos = halos[list('xyz')].view((float, 3))

for j, m_log in zip(k, m_log_list):
    fname = 'HOD_pair_counts/{}_nd{:.1f}'.format(case, m_log)
    np.save(fname, count_pairs(halos[s[j:]], rbins, float(box_size)))
