#!/usr/bin/env python

from itertools import product
from subprocess import check_call

proxies = ('vpeak', 'vmax', 'macc')
cases = ('c250-2048', 'c250-1024', 'c250-768', 'c250-512', 'c400-1024', 'c400-768', 'c125-1024')
#proxies = ('mvir',)
#cases = ('c250-2048', 'c250-1024', 'c250-768', 'c250-512', 'c400-2048', 'c400-1024', 'c400-768', 'c125-2048', 'c125-1024')
#cases = ('c250-2560',)
run = 0

for proxy, case in product(proxies, cases):
#     check_call(['bsub.new', '-m', 'bulletfarm', '-W', '40:00', '-R', '"mem > 30000"','-M', '64000', '~/bin/python-wait', 'pairs.py', proxy, case])
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '10:00', '~/bin/python-wait', 'save_simulate_lowres_mass_function.py', proxy, case, '-r', str(run)])