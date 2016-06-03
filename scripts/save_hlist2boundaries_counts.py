import numpy as np
from helpers.SimulationAnalysis import readHlist
from fast3tree import fast3tree
from argparse import ArgumentParser
import os.path

def point2JKindex(loc_array, box_size = 125., n = 50):
    dL = box_size/float(n)
    index = np.floor_divide(loc_array, dL).astype(int)
    return np.ravel_multi_index(index.T, (n,n,n))
    
def JKindex2JKbox(idx, box_size = 125., n = 50):
    i,j,k = np.unravel_index(idx, (n,n,n))
    dL = float(box_size)/n
    return dL*np.array([[i, i+1],[j, j+1], [k, k+1]])

def JKindex2boxcenter(idx, box_size = 125., n = 50):
    i,j,k = np.unravel_index(idx, (n,n,n))
    dL = float(box_size)/n
    return dL*np.array([i+0.5, j+0.5, k+0.5])


parser = ArgumentParser()
parser.add_argument('proxy')
parser.add_argument('case')
parser.add_argument('--nsave', type = int, default = 12)
args = parser.parse_args()
proxy = args.proxy
case = args.case
n = args.nsave

#settings
nd_log_list = np.linspace(-3.3, -1.7, 17)
#nd_log_list = np.linspace(-2.3, -1.8, 6)
#nd_log_list = np.linspace(-1.7, -3.3, 17)
rbins = np.logspace(-1.3, 1.3, 27)
box_size = float(case[1:4])

#import data
hlist_path ='/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case)
halos_all = np.load(hlist_path)
s = halos_all[proxy].argsort()
kk = (10.0**nd_log_list * (box_size**3)).astype(int)
kk *= -1
halos_all = halos_all[list('xyz')].view((float, 3))


for jj, nd_log in zip(kk, nd_log_list):
    # select the halos above the nd limit
    halos = halos_all[s[jj:]]
    rbins = np.logspace(-1.3, 1.3, 27)
    #JKcrosspairs will save the number of halo partners that are inside each of the 27 neighboring boxes
    JKcrosspairs = np.zeros((n**3, len(rbins), 27))
    for JKidx in range(n**3):
        # select the points inside the JK box
        xr, yr, zr = JKindex2JKbox(JKidx, box_size, n)
        xbool = (xr[0] - halos[:,0])*(xr[1] - halos[:,0]) < 0
        ybool = (yr[0] - halos[:,1])*(yr[1] - halos[:,1]) < 0
        zbool = (zr[0] - halos[:,2])*(zr[1] - halos[:,2]) < 0
        halos_JK = halos[np.all([xbool, ybool, zbool],axis=0)]

        #determine the neighboring 27 JK indices (including itself)
        center = np.mean([xr, yr, zr], axis=-1)
        xx, yy, zz = [grid.flatten() for grid in np.meshgrid([-1.,0,1.],[-1.,0,1.],[-1.,0,1.],indexing='ij')]
        dx, dy, dz = float(box_size)/n*np.identity(3)
        shifts = np.outer(xx, dx) + np.outer(yy, dy) + np.outer(zz,dz)
        neighbor_JKindex = point2JKindex(np.mod(center + shifts , box_size), box_size, n)
        # calculate the number of halo partners that are inside each of the 27 neighboring boxes
        # crosspairslist[r, 13] is the number of autopairs
        # start the autopairs with negative so that we don't count the halo pairing with itself
        crosspairslist = np.zeros((len(rbins), 27))
        crosspairslist[:,13] = -len(halos_JK)*np.ones_like(rbins)
        with fast3tree(halos) as tree:
            tree.set_boundaries(0, box_size)
            for p in halos_JK:
                for i, r in enumerate(rbins):
                    pos = tree.query_radius(p,r, periodic=True, output='pos')
                    JKpos = point2JKindex(pos, box_size, n)
                    for k, neighborJK in enumerate(neighbor_JKindex):
                        crosspairslist[i,k] += np.count_nonzero(JKpos == neighborJK)
        JKcrosspairs[JKidx] = crosspairslist
    fname = 'pairs_{}/boundaries/{}_nd{:.1f}_nsave{}'.format(proxy, case, nd_log, n)
    np.save(fname, JKcrosspairs)

