#!/usr/bin/env python

from itertools import product
from subprocess import check_call

proxies = ('vpeak', 'vmax', 'mvir','macc')
proxies = ('vpeak','vmax','macc')
case = 'MDhigh'
nside = 8
delta = 64
idxs = range(0,512,delta)

for proxy, idx in product(proxies, idxs):
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '40:00', '~/bin/python-wait', 'save_remove_hosthalos.py', proxy, str(idx), str(idx + delta), '-n', str(nside), '-c', case])