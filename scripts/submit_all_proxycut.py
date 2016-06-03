#!/usr/bin/env python

from itertools import product
from subprocess import check_call

proxies = ('vpeak', 'vmax', 'mvir')
case = 'darksky'
nside = 4
delta = 8
idxs = range(0,64,delta)

for proxy, idx in product(proxies, idxs):
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '02:00', '~/bin/python-wait', 'save_proxycut.py', proxy, str(idx), str(idx + delta), '--nside', str(nside), '-c', case])