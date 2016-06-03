#!/usr/bin/env python

from itertools import product
from subprocess import check_call

proxies = ('vpeak', 'vmax', 'mvir', 'macc')
#proxies = ('vpeak', 'vmax', 'macc')
nside = 9

delta = 81
#idxs = [18,9,0]
idxs = range(0,729,delta)

for proxy, idx in product(proxies, idxs):
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '20:00', '~/bin/python-wait', 'save_pair_counts_multidark_all.py', proxy, str(idx), str(idx + delta), '-n', str(nside)])
    
#    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '20:00', '~/bin/python-wait', 'save_pair_counts_multidark_all.py', proxy, str(idx), str(idx + delta), '--nside', str(nside), '-b', '1000', '-p', '/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study/multidark/hlist_0.68220.npy'])
    
#    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '40:00', '~/bin/python-wait', 'save_pair_counts_multidark_all.py', proxy, str(idx), str(idx + delta), '-n', str(nside), '-b', '400' , '-p', '/u/ki/yymao/ki-des/ds14_i_4096/hlist_1.00000.npy'])

#    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '40:00', '~/bin/python-wait', 'save_pair_counts_multidark_all.py', proxy, str(idx), str(idx + delta), '-n', str(nside), '-b', '400' , '-p', '/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/c400-2048/hlist_1.00000.npy'])