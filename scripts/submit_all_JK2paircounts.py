#!/usr/bin/env python

from itertools import product
from subprocess import check_call

proxies = ('vpeak', 'vmax', 'mvir', 'macc')
#cases = ('c125-2048-new',)
#proxies = ('vmax',)
#cases = ('c250-2048', 'c250-1024', 'c250-768', 'c250-512', 'c400-2048', 'c400-1024', 'c400-768', 'c125-2048', 'c125-1024')
cases = ('c400-2048', 'c400-1024', 'c400-768')
#proxies = ('vmax', 'mvir')
#cases = ('rock-1', 'rock-2', 'rock-3.5', 'rock-5')
nJK = str(5)

for proxy, case in product(proxies, cases):
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '01:00', '~/bin/python-wait', 'JKindex2paircounts.py', proxy, case, '--nJK', nJK])

