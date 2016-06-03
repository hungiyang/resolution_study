__all__=['af', 'urlopen', 'AbundanceFunction', 'filepath', 'os', 'np', '_projected_samples_cache', 'load_projected_samples', 'load_precalculated_JK', '_pairs_cache_JK', 'load_sim_pair_count', '_pairs_sim_cache', 'load_samples', '_samples_cache', 'load_pair_count', '_pairs_cache', 'load_pair_count_HOD', '_pairs_HOD_cache', 'load_nd_limits', 'load_proxy_cut', '_proxycut_cache', '_nd_limits', 'parse_boxname', 'frexp10']

import filepath
import os
import numpy as np
from AbundanceMatching import AbundanceFunction
from urllib import urlopen

af = np.loadtxt(urlopen('https://arxiv.org/src/1304.7778v2/anc/LF_Ser.dat'), usecols=(0,1))
af[:,0] -= 5.0*np.log10(0.7)
af[:,1] = (10.**af[:,1])/(0.7**3)
af = AbundanceFunction(*af.T, ext_range=(-25, -16), faint_end_fit_points=5)


_projected_samples_cache={}
def load_projected_samples(proxy, nd_log, box_size = 125., case = 'MDhigh', verbal = False):
    key = (proxy, nd_log, box_size, case)
    if key not in _projected_samples_cache:
        if case == 'MDhigh':
            _projected_samples_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_highres + 'projected/boxsize_{}/pairs_{}/nd{:.1f}_idx{:03d}.npy'.format(int(box_size),proxy,nd_log, i)
                if os.path.isfile(fname):
                    _projected_samples_cache[key].extend(np.load(fname).reshape((1,-1)))
            _projected_samples_cache[key] = np.array(_projected_samples_cache[key])
            if verbal: print 'loaded{}{}{}{}'.format(*key)
            return _projected_samples_cache[key]
    return _projected_samples_cache[key]


_pairs_cache_JK = {}
def load_precalculated_JK(proxy, case, nd_log, nJK = 5):
    key = (proxy, nJK, case, nd_log)
    if key not in _pairs_cache_JK:
        _pairs_cache_JK[key] = np.load(filepath.path_jackknife + 'pairs_{}/nJK_{}/{}_nd{:.1f}.npy'.format(*key))    
    return np.sum(_pairs_cache_JK[key].reshape((nJK**3,-1,2)), axis = -1)
    #return _pairs_cache_JK[key]

_pairs_sim_cache = {}
def load_sim_pair_count(proxy, case, nd_log, cumulative=False):
    key = (proxy, case, nd_log)
    if key not in _pairs_cache:
        if os.path.isfile(filepath.path_sim_lowres_counts + '/pairs_{}/{}_nd{:.1f}_run0.npy'.format(*key)):
            _pairs_sim_cache[key] = np.load(filepath.path_sim_lowres_counts + 'pairs_{}/{}_nd{:.1f}_run0.npy'.format(*key)).astype(float)      
        else:
            _pairs_sim_cache[key] = np.load(filepath.path_sim_lowres_counts + 'pairs_{}/{}_nd{:.1f}_run1.npy'.format(*key)).astype(float)      
    return _pairs_sim_cache[key] if cumulative else np.ediff1d(_pairs_sim_cache[key][::2])
    
    
_samples_cache={}
def load_samples(proxy, nd_log, box_size = 125., case = 'MDhigh', verbal = False):
    key = (proxy, nd_log, box_size, case)
    if key not in _samples_cache:
        if case == 'darksky':
            _samples_cache[key] = []
            n = int(400/box_size)
            for i in range(n**3):
                fname = filepath.path_dark_sky+'pairs_{}/boxsize_{}/nd{:.1f}_idx{:03d}.npy'.format(proxy,int(box_size) ,nd_log, i)
                if os.path.isfile(fname):
                    _samples_cache[key].extend(np.ediff1d(np.load(fname)).reshape((1,-1)))
            _samples_cache[key] = np.array(_samples_cache[key])
            if verbal: print 'loaded{}{}{}{}'.format(*key)
            return _samples_cache[key]
        
        if case == 'MDlow':
            _samples_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_lowres + 'pairs_{}/boxsize_{}/nd{:.1f}_idx{:03d}.npy'.format(proxy,int(box_size),nd_log , i)
                if os.path.isfile(fname):
                    _samples_cache[key].extend(np.ediff1d(np.load(fname)).reshape((1,-1)))
            _samples_cache[key] = np.array(_samples_cache[key])
            if verbal: print 'loaded{}{}{}{}'.format(*key)
            return _samples_cache[key]
        # save the combined resultes directly to make loading faster
        if case == 'MDhigh':
            fname = filepath.path_multidark + 'load_samples_combined/pairs_{}/bs{:d}_nd{:.1f}.npy'.format(proxy, int(box_size), nd_log)
            _samples_cache[key] = np.load(fname)
            return _samples_cache[key]
        '''
        if case == 'MDhigh':
            _samples_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_highres + 'pairs_{}/boxsize_{}/nd{:.1f}_idx{:03d}.npy'.format(proxy,int(box_size),nd_log , i)
                if os.path.isfile(fname):
                    _samples_cache[key].extend(np.ediff1d(np.load(fname)).reshape((1,-1)))
            _samples_cache[key] = np.array(_samples_cache[key])
            if verbal: print 'loaded{}{}{}{}'.format(*key)
            return _samples_cache[key]
        '''
        
        if case == 'c400-2048' or case=='lowres_z1' or case=='lowres_z3':
            _samples_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark + '{}/pairs_{}/boxsize_{}/nd{:.1f}_idx{:03d}.npy'.format(case, proxy,int(box_size),nd_log , i)
                if os.path.isfile(fname):
                    _samples_cache[key].extend(np.ediff1d(np.load(fname)).reshape((1,-1)))
            _samples_cache[key] = np.array(_samples_cache[key])
            if verbal: print 'loaded{}{}{}{}'.format(*key)
            return _samples_cache[key]
            
    return _samples_cache[key]

_pairs_cache = {}
def load_pair_count(proxy, case, nd_log, cumulative=False, fixed_limits=False):
    key = (proxy, case, nd_log, 'f' if fixed_limits else 'c')
    if key not in _pairs_cache:
        _pairs_cache[key] = np.load(filepath.path_simple_pair_counts + 'pairs_{}/{}_nd{:.1f}_{}.npy'.format(*key)).astype(float)      
    return _pairs_cache[key] if cumulative else np.ediff1d(_pairs_cache[key][::2])

_pairs_HOD_cache = {}
def load_pair_count_HOD(case, m_log, cumulative=False):
    key = (case, m_log)
    if key not in _pairs_cache:
        _pairs_HOD_cache[key] = np.load(filepath.path_HOD_pair_counts + '{}_nd{:.1f}.npy'.format(*key)).astype(float)      
    return _pairs_HOD_cache[key] if cumulative else np.ediff1d(_pairs_HOD_cache[key][::2])

_nd_limits = {}
def load_nd_limits(proxy, case):
    key = (proxy, case)
    if key not in _nd_limits:
        _nd_limits[key] = np.load(filepath.path_nd_limits + 'pairs_{}/{}_nd_limits.npy'.format(*key))     
    return _nd_limits[key]

_proxycut_cache = {}
def load_proxy_cut(proxy, box_size, case):
    key = (proxy, box_size, case)
    if key not in _proxycut_cache:
        if case == 'MDhigh':
            _proxycut_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_highres + 'pairs_{}/boxsize_{}/proxycut_idx{:03d}.npy'.format(proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        
        #extended nd_log_list
        if case == 'MDhigh_more':
            _proxycut_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_highres + 'pairs_{}/boxsize_{}/proxycut/proxycut_idx{:03d}.npy'.format(proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
    
        if case == 'MDlow':
            _proxycut_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark_lowres + 'pairs_{}/boxsize_{}/proxycut_idx{:03d}.npy'.format(proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        if case == 'lowres_z1' or case == 'lowres_z3':
            _proxycut_cache[key] = []
            n = int(1000/box_size)
            for i in range(n**3):
                fname = filepath.path_multidark + '{}/pairs_{}/boxsize_{}/proxycut_idx{:03d}.npy'.format(case,proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        
        if case == 'darksky':
            _proxycut_cache[key] = []
            n = int(400/box_size)
            for i in range(n**3):
                fname = filepath.path_dark_sky + 'pairs_{}/boxsize_{}/proxycut_idx{:03d}.npy'.format(proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        
        if case == 'darksky_more':
            _proxycut_cache[key] = []
            n = int(400/box_size)
            for i in range(n**3):
                fname = filepath.path_dark_sky + 'pairs_{}/boxsize_{}/proxycut/proxycut_idx{:03d}.npy'.format(proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        
        if case == 'c400-2048':
            _proxycut_cache[key] = []
            n = int(400/box_size)
            for i in range(n**3):
                fname = filepath.path_dark_sky + '{}/pairs_{}/boxsize_{}/proxycut_idx{:03d}.npy'.format(case, proxy,int(box_size), i)
                if os.path.isfile(fname):
                    _proxycut_cache[key].extend(np.load(fname).reshape((1,-1)))
            _proxycut_cache[key] = np.array(_proxycut_cache[key])
            return _proxycut_cache[key]
        
        #_proxycut_cache[proxy] = np.zeros((512,17), float)
        #for i in range(512):
        #    _proxycut_cache[proxy][i] = np.load('/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/pairs_mvir/boxsize_{}/proxycut_idx{:03d}.npy'.format(proxy, i))
    return _proxycut_cache[key]

def parse_boxname(boxname):
    items = boxname[1:].split('-')
    box_size = float(items[0])
    npart = int(items[1])
    return box_size, npart, (box_size/npart)**3*79379292719.00159

def frexp10(x):
    e = int('{0:e}'.format(float(x)).partition('e')[-1])
    return float(x)/(10.0**e), e
