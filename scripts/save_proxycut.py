from argparse import ArgumentParser
import numpy as np
from fast3tree import fast3tree

parser = ArgumentParser()
parser.add_argument('proxy')
parser.add_argument('idxmin', type = int, default = 0)
parser.add_argument('idxmax', type = int, default = 3**3)
parser.add_argument('-n','--nside', type = int, default = 8)
parser.add_argument('-c', '--case', default = 'MDhigh')
args = parser.parse_args()
proxy = args.proxy
idxmin = args.idxmin
idxmax = args.idxmax
box_size_all = 1000.

nd_log_list = np.linspace(-3.3, -2.0, 17)
case = args.case
if case=='MDhigh': 
    hlist_path = '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy'
    nd_log_list = np.linspace(-3.3, -1.1, 23)
if case=='darksky': 
    hlist_path = '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy'
    box_size_all = 400.
    nd_log_list = np.linspace(-3.3, -0.5, 29)
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

#hlist_path = '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy'
halos_all = np.load(hlist_path)[['x','y','z',proxy]]

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

    points = halos[list('xyz')].view((float, 3))

    #save the threshold proxy value at each nd_log
    proxy_cut = np.fromiter((halos[proxy][s[j]] for j in k), float)
    #high res Multidark
    if hlist_path == '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy':
        fname = 'multidark/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    #low res Multidark
    if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_1.00110.npy':
        fname = 'multidark/lowres/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    #low res Multidark z=1
    if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.49990.npy':
        fname = 'multidark/lowres_z1/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    #low res Multidark z=3
    if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.25690.npy':
        fname = 'multidark/lowres_z3/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    #darksky 400-4096 box
    if hlist_path == '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy':
        fname = 'darksky/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    #c400-2048 box
    if hlist_path == '/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/c400-2048/hlist_1.00000.npy':
        fname = 'multidark/c400-2048/pairs_{}/boxsize_{:d}/proxycut/proxycut_idx{:03d}'.format(proxy, int(box_size), idx)
    np.save(fname, proxy_cut)
    

