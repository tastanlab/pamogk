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
            print(f'  Writing kernel file {i + 1:3}/{s[0]}', end='\r')
            with open(f'/tmp/kernel-{i}.csv', 'w') as f:
                sio.savemat(f'/tmp/kms-{i}', {f'out_kernel{i}_': m})
                c = csv.writer(f)
                c.writerows(m)
        sio.savemat('/tmp/kms.mat', {'kms': kms})
        print('\nFinished')
    except:
        print('Failed to convert kernel matrices from given path:', p)
        traceback.print_exc()
