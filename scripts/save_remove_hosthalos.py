from argparse import ArgumentParser
import numpy as np
from fast3tree import fast3tree

parser = ArgumentParser()
parser.add_argument('proxy')
parser.add_argument('idxmin', type = int, default = 0)
parser.add_argument('idxmax', type = int, default = 3**3)
parser.add_argument('-n','--nside', type = int, default = 8)
#right now lowres Multidark cases don't work because they don't have id, rvir fields.
parser.add_argument('-c', '--case', default = 'MDhigh')
parser.add_argument('-m','--mcut', type = float, default = 1.e14)
args = parser.parse_args()
proxy = args.proxy
idxmin = args.idxmin
idxmax = args.idxmax
mcut = args.mcut
box_size_all = 1000.

#settings
rbins = np.logspace(-1.3, 1.3, 27)
nd_log_list = np.linspace(-3.3, -2.0, 17)
case = args.case
if case=='MDhigh': 
    hlist_path = '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy'
    nd_log_list = np.linspace(-3.3, -1.7, 17)
if case=='darksky': 
    hlist_path = '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy'
    box_size_all = 400.
    nd_log_list = np.linspace(-3.3, -1.7, 17)
if case=='MDlow_z1': 
    hlist_path = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.49990.npy'
if case=='MDlow_z3': 
    hlist_path = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.25690.npy'
if case=='MDlow': 
    hlist_path = '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy'

box_size = box_size_all/args.nside

#select_sample takes the halo catalog and return the halos inside the smaller sample box
def select_sample(halos, idx, box_size):
    pos = np.mod(halos[['x','y','z']].view((float, 3)),box_size_all)
    n = int(box_size_all/box_size)
    halos_idx = halos[np.ravel_multi_index(np.floor_divide(pos,box_size).astype(int).T, (n,n,n)) == idx]
    return halos_idx

# count number of pairs in each rbin radius (cumulative)
def count_pairs(points, rbins):
    pairs = np.zeros(len(rbins), dtype=int)
    with fast3tree(points) as tree:
        #tree.set_boundaries(0, box_size)
        #no periodic boundary condition
        for p in points:
            pairs += np.array([tree.query_radius(p, r, periodic=False, \
                    output='c') for r in rbins])
    return (pairs-points.shape[0])/2

def remove_subhalos(halos, halos_nd, mcut = 1.e14):
    hosts = halos[halos['mvir']>mcut][['x','y','z']].view((float,3))
    rvir = halos[halos['mvir']>mcut][['rvir']].view((float,1))
    hostidx = halos[halos['mvir']>mcut][['id']].view((int,1))
    hostidx_nd = np.array([np.where(halos_nd['id'].view((int,1))==idx)[0] for idx in hostidx])
    #deal with the case that the host halos are not in the selected nd sample
    hostidx_nd = np.array([idx[0] for idx in hostidx_nd if len(idx)==1])
    points = halos_nd[list('xyz')].view((float,3))
    
    not_subhalo_bool = np.ones(len(points), dtype=int)
    with fast3tree(points) as tree:
        #tree.set_boundaries(0, box_size)
        #no periodic boundary condition
        for i, (p, rv, idx_nd) in enumerate(zip(hosts, rvir, hostidx_nd)):
            # rv is in unit Kpc, so divide by 1000 to get units in Mpc
            sublist = tree.query_radius(p, 2*rv/1000., periodic=False, output='index')
            for k in sublist:
                if k!= idx_nd: 
                    not_subhalo_bool[k] = 0
    return halos_nd[np.where(not_subhalo_bool)]

if proxy != 'mvir':
    halos_all = np.load(hlist_path)[['x','y','z',proxy,'id','rvir','mvir']]
else:
    halos_all = np.load(hlist_path)[['x','y','z',proxy,'id','rvir']]

# add the option of going from large to small index
if idxmin > idxmax:
    idxrange = range(idxmax, idxmin)[::-1]
else: 
    idxrange = range(idxmin, idxmax)
    
for idx in idxrange:
    #run
    halos = select_sample(halos_all, idx, box_size)
    s = halos[proxy].argsort()

    k = (10.0**nd_log_list * (float(box_size)**3)).astype(int)
    k *= -1

    for j, nd_log in zip(k, nd_log_list):
        points = remove_subhalos(halos, halos[s[j:]], mcut)[['x','y','z']].view((float,3))
        
        #high res Multidark
        if hlist_path == '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy':
            fname = 'multidark/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        #low res Multidark
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_1.00110.npy':
            fname = 'multidark/lowres/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        #low res Multidark z=1
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.49990.npy':
            fname = 'multidark/lowres_z1/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        #low res Multidark z=3
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.25690.npy':
            fname = 'multidark/lowres_z3/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        #darksky 400-4096 box
        if hlist_path == '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy':
            fname = 'darksky/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        #c400-2048 box
        if hlist_path == '/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/c400-2048/hlist_1.00000.npy':
            fname = 'multidark/c400-2048/pairs_{}/boxsize_{}/remove_host/nd{:.1f}_idx{:03d}'.format(proxy, int(box_size), nd_log, idx)
        np.save(fname, count_pairs(points, rbins))
    