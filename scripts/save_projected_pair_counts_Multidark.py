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


# from Yao's helpers.CorrelationFunction
def _yield_periodic_points(center, dcorner1, dcorner2, box_size):
    cc = np.array(center)
    flag = (cc+dcorner1 < 0).astype(int) - (cc+dcorner2 >= box_size).astype(int)
    cp = cc + flag*box_size
    a = range(len(cc))
    for j in xrange(1 << len(cc)):
        for i in a:
            if j >> i & 1 == 0:
                cc[i] = center[i]
            elif flag[i]:
                cc[i] = cp[i]
            else:
                break
        else:
            yield cc

# code of projected correlation function from Yao's helper module
def projected_correlation(points, rbins, zmax , box_size):
    points = np.asarray(points)
    s = points.shape
    if len(s) != 2 or s[1] != 3:
        raise ValueError('`points` must be a 2-d array with last dim=3')
    N = s[0]

    rbins = np.asarray(rbins)
    rbins_sq = rbins*rbins
    dcorner2 = np.array([rbins[-1], rbins[-1], zmax])
    dcorner1 = -dcorner2
    if np.any(dcorner2*2 > box_size):
        print "[Warning] box too small!"

    pairs_rand = float(N*N) / box_size**3 \
            * (rbins[1:]**2-rbins[:-1]**2)*np.pi*zmax*2.0

    dcorner1[2] = 0 #save some time
    pairs = np.zeros(len(rbins)-1, dtype=int)
    with fast3tree(points) as tree:
        for p in points:
            for pp in _yield_periodic_points(p,dcorner1,dcorner2,box_size):
                x,y=tree.query_box(pp+dcorner1,pp+dcorner2,output='p').T[:2]
                x -= pp[0]; x *= x
                y -= pp[1]; y *= y
                x += y; x.sort()
                pairs += np.ediff1d(np.searchsorted(x, rbins_sq))
    return pairs

#select_sample takes the halo catalog and return the halos inside the smaller sample box
def select_sample(halos, idx, box_size):
    pos = np.mod(halos[['x','y','z']].view((float, 3)),box_size_all)
    n = int(box_size_all/box_size)
    halos_idx = halos[np.ravel_multi_index(np.floor_divide(pos,box_size).astype(int).T, (n,n,n)) == idx]
    return halos_idx

#settings
box_size_all = 1000.
nd_log_list = np.linspace(-3.3, -1.7, 17)
rbins = np.logspace(-1.3, 1.3, 27)
case = args.case
if case=='MDhigh': 
    hlist_path = '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy'
if case=='darksky': 
    hlist_path = '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy'
    box_size_all = 400.
if case=='MDlow_z1': 
    hlist_path = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.49990.npy'
    nd_log_list = np.linspace(-3.3, -2.0, 14)
if case=='MDlow_z3': 
    hlist_path = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.25690.npy'
    nd_log_list = np.linspace(-3.3, -2.2, 12)
if case=='MDlow': 
    hlist_path = '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_1.00110.npy'
    nd_log_list = np.linspace(-3.3, -2.0, 14)
box_size = box_size_all/args.nside

# import data
halos_all = np.load(hlist_path)[['x','y','z','vz',proxy]]

# add the option of going from large to small index
if idxmin > idxmax:
    idxrange = range(idxmax, idxmin)[::-1]
else: 
    idxrange = range(idxmin, idxmax)
    
for idx in idxrange:
    #run
    halos = select_sample(halos_all, idx, box_size)
    #add redshift distortion
    halos['z'] = halos['z'] + halos['vz']/100.
    s = halos[proxy].argsort()

    k = (10.0**nd_log_list * (float(box_size)**3)).astype(int)
    k *= -1

    points = halos[list('xyz')].view((float, 3))

    for j, nd_log in zip(k, nd_log_list):
        #high res Multidark
        if hlist_path == '/u/ki/lehmann/lustre/multidark/NewMDPL/hlists/hlist_1.00000.npy':
            fname = 'multidark/highres/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        #low res Multidark
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_1.00110.npy':
            fname = 'multidark/lowres/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        #low res Multidark z=1
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.49990.npy':
            fname = 'multidark/lowres_z1/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        #low res Multidark z=3
        if hlist_path == '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.25690.npy':
            fname = 'multidark/lowres_z3/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        #darksky 400-4096 box
        if hlist_path == '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy':
            fname = 'darksky/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        #c400-2048 box
        if hlist_path == '/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/c400-2048/hlist_1.00000.npy':
            fname = 'multidark/c400-2048/projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}'.format(int(box_size), proxy,  nd_log, idx)
        np.save(fname, projected_correlation(points[s[j:]], rbins, 40., box_size))
