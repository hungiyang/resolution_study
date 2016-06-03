from argparse import ArgumentParser
import numpy as np
from fast3tree import fast3tree

parser = ArgumentParser()
parser.add_argument('proxy')
parser.add_argument('case')
parser.add_argument('-r','--run', type= int, default = 0)
args = parser.parse_args()
case = args.case
proxy = args.proxy
run = args.run

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
#nd_log_list = np.linspace(-3.3, -1.7, 17)
nd_log_list = np.linspace(-5.0,-1.7,34)

box_size = case[1:4]
hlist_path_low ='/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case)
if float(box_size) == 250.:
    case_high = case.partition('-')[0] + '-2560'
else:
    case_high = case.partition('-')[0] + '-2048'
hlist_path_high ='/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case_high)

#run
halos_high = np.load(hlist_path_high)
halos_low = np.load(hlist_path_high)

#convert nd limits to the proxy limits of the low res box
k = (10.**nd_log_list*float(box_size)**3).astype(int)
s = halos_low[proxy].argsort()
proxy_log_list = np.array([np.log10(halos_low[proxy][s[-j]]) for j in k])
proxy_log_list = np.insert(proxy_log_list,0, np.inf)

# calculate the number of halos in each proxy bin in the lowres case
n_halos_list = np.insert(np.ediff1d(k) ,0, k[0])


# select random sample from high resolution
sim_idx = [[] for i in range(len(n_halos_list))]
temp = np.array([]).astype(int)
for i,(high, low) in enumerate(zip(proxy_log_list[:-1], proxy_log_list[1:])):
    temp = np.concatenate((temp, np.random.choice(np.where(np.logical_and(halos_high[proxy] > 10.**low,halos_high[proxy] < 10.**high))[0], n_halos_list[i])))
    sim_idx[i] = temp
sim_idx = np.array(sim_idx)

halos = halos_high[list('xyz')].view((float, 3))

ndidxcut = np.abs(nd_log_list - -3.3).argmin()
for nd_log, idxlist in zip(nd_log_list, sim_idx)[ndidxcut:]:
    fname = 'random_sample_pair_counts/pairs_{}/{}_nd{:.1f}_run{}'.format(proxy, case, nd_log, run)
    np.save(fname, count_pairs(halos[idxlist], rbins, float(box_size)))
