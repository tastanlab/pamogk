#!/usr/bin/env python3
import csv
import sys
import traceback

import numpy as np
import scipy.io as sio

args = sys.argv[1:]

for p in args:
    try:
        print('Converting kernel matrices from file:', p)
        kms = np.load(p)['kms']
        s = kms.shape
        i_lower = np.tril_indices(s[1], -1)
        for i, m in enumerate(kms):
            m[i_lower] = m.T[i_lower]  # for precision errors mirror upper half
            print('  Writing kernel file {:3}/{}'.format(i + 1, s[0]), end='\r')
            with open('/tmp/kernel-{}.csv'.format(i), 'w') as f:
                sio.savemat('/tmp/kms-{}'.format(i), {'out_kernel{}_'.format(i): m})
                c = csv.writer(f)
                for r in m: c.writerow(r)
        sio.savemat('/tmp/kms.mat', {'kms': kms})
        print('\nFinished')
    except:
        print('Failed to convert kernel matrices from given path:', p)
        traceback.print_exc()
