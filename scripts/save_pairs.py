from argparse import ArgumentParser
import numpy as np
from count_pairs_JKindex import count_pairs_JKindex
from count_pairs import count_pairs

parser = ArgumentParser()
parser.add_argument('proxy')
parser.add_argument('case')
parser.add_argument('--JK', action='store_true')
parser.add_argument('--fixed-limits', action='store_true')
args = parser.parse_args()
proxy = args.proxy
case = args.case

#settings
rbins = np.logspace(-1.3, 1.3, 27)
#nd_log_list = np.linspace(-3.3, -1.7, 17)
nd_log_list = np.linspace(-1.7, -3.3, 17)
#nd_log_list = np.linspace(-2.0, -1.7, 4)
out_dir = 'pairs_' + proxy

if case.startswith('rock-'):
    box_size = 125.0
    hlist_path ='/u/ki/yymao/ki21/sham_test/resolution-test/c125-1024/rockstar/out_99_s{}.npy'.format(case[5:])
else:
    box_size = case[1:4]
    hlist_path ='/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case)

#run
halos = np.load(hlist_path)
s = halos[proxy].argsort()

if args.fixed_limits:
    limits = np.load('pairs_{}/c{}_nd_limits.npy'.format(proxy, box_size))
    k = np.searchsorted(halos[proxy], limits, sorter=s)
else:
    k = (10.0**nd_log_list * (float(box_size)**3)).astype(int)
    k *= -1

halos = halos[list('xyz')].view((float, 3))

if args.JK:
    for j, nd_log in zip(k, nd_log_list):
        fname = 'pairs_{}/{}_nd{:.1f}_{}'.format(proxy, case, nd_log, \
                'f' if args.fixed_limits else 'c')
        np.save(fname, count_pairs_JKindex(halos[s[j:]], rbins, float(box_size)))
        
else:
    for j, nd_log in zip(k, nd_log_list):
        fname = 'simple_pair_counts/pairs_{}/{}_nd{:.1f}_{}'.format(proxy, case, nd_log, \
                'f' if args.fixed_limits else 'c')
        np.save(fname, count_pairs(halos[s[j:]], rbins, float(box_size)))

