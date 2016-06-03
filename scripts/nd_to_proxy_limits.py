import numpy as np

proxies = ('vpeak', 'vmax', 'mvir', 'macc')
box_sizes = (125, 250, 400)
cases = ('c250-2048', 'c250-1024', 'c250-768', 'c250-512', 'c400-2048', 'c400-1024', 'c400-768', 'c125-2048', 'c125-1024')
nd_log_list = np.linspace(-3.3, -1.7, 17)

for case in cases:
    halos = np.load('/u/ki/yymao/ki21/sham_test/resolution-test/{}/hlist_1.00000.npy'.format(case))
    k = np.around(10.0**nd_log_list * (float(case[1:4])**3)).astype(int)
    flag = (k > len(halos))
    k[flag] = 0
    k[~flag] *= -1
    for proxy in proxies:
        s = halos[proxy].argsort()
        np.save('pairs_{}/{}_nd_limits'.format(proxy, case), \
                halos[proxy][s[k]])

