import numpy as np
from argparse import ArgumentParser
import os.path

def pair_count_JK(pairs_list, JKindex, nJK = 2, nsaved = 50):
    pairs_JK = np.zeros(len(pairs_list),'float')
    idx, idy, idz = np.unravel_index(np.arange(nsaved**3),(nsaved, nsaved, nsaved))
    r = int(nsaved/nJK)
    index_conversion = np.ravel_multi_index((idx/r, idy/r, idz/r),(nJK, nJK, nJK))
    index2bool = index_conversion != JKindex
    for i, pairs in enumerate(pairs_list):
        # calculate the number of pairs that are not in the JK box (both indices outside of the JK box)
        if len(pairs) == 0: pairs_JK[i] = 0
        else:
            JKbool = index2bool[pairs]
            pairs_JK[i] = np.count_nonzero(np.all(JKbool, axis = 1))
    return np.ediff1d(pairs_JK)

_pairs_cache_JK = {}
def load_pair_countJK(proxy, case, nd_log, fixed_limits=False):
    key = (proxy, case, nd_log, 'f' if fixed_limits else 'c')
    if key not in _pairs_cache_JK:
        fn = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/pairs_{}/{}_nd{:.1f}_{}.npy'.format(*key)
        _pairs_cache_JK[key] = np.load(fn)    
    return _pairs_cache_JK[key]

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('proxy')
    parser.add_argument('case')
    parser.add_argument('--nJK', type = int, default = 5)
    parser.add_argument('--nsaved', type = int, default = 50)
    args = parser.parse_args()
    proxy = args.proxy
    case = args.case
    nJK = args.nJK
    nsaved = args.nsaved
    
    nd_log_list = np.linspace(-3.3, -1.7, 17)
#     nd_log_list = np.linspace(-2.3, -2.3, 1)
    
    for nd_log in nd_log_list:
        key = (proxy, case, nd_log, 'c')
        fn = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/pairs_{}/{}_nd{:.1f}_{}.npy'.format(*key)
        if os.path.isfile(fn):
            pairJKindex = load_pair_countJK(proxy, case, nd_log)
            paircounts = np.array([pair_count_JK(pairJKindex, i, nJK = nJK, nsaved = nsaved) for i in range(nJK**3)])
            fname = 'pairs_{}/nJK_{}/{}_nd{:.1f}'.format(proxy, nJK, case, nd_log)
            np.save(fname, paircounts)






