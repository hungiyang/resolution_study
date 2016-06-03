#!/usr/bin/env python
from itertools import product
from subprocess import check_call

cases = ('c250-2560','c250-2048', 'c250-1024', 'c250-768', 'c250-512', 'c400-2048', 'c400-1024', 'c400-768', 'c125-2048', 'c125-1024')

for case in cases:
    check_call(['bsub.new', '-m', 'bulletfarm', '-W', '05:00', '~/bin/python-wait', 'save_pair_counts_HOD.py', case])