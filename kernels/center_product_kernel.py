import numpy as np
from sklearn.svm import SVC
import math

def calculate_S_and_P(patients, gene_vectors):
    # calculate S (mutated gene vector set) and P (average mutataion point) vector
    for p in patients:
        genes = np.array([gene_vectors[str(n)] for n in p['mutated_nodes']])
        p['S'] = genes
        p['P'] = np.average(genes, axis=0)

def test_accr(patients):
    hit = 0
    pids = [p['pid'] for p in patients]
    for pid in pids:
        test_p = [p for p in patients if p['pid'] == pid][0]
        train_p = [p for p in patients if p['pid'] != pid]

        linear_svc = SVC(kernel='linear')
        linear_svc.fit([p['P'] for p in train_p], [p['sick'] for p in train_p])
        is_hit = linear_svc.predict([test_p['P']]) == [test_p['sick']]
        # print('%3d: %s' % (pid, is_hit), linear_svc.predict([test_p['P']]), test_p['pid'], test_p['sick'])
        hit += is_hit[0]
    print('Accuracy Leave-One-Out accuracy=%.2lf' % (hit/len(pids)))

    hit = 0
    pids = [p['pid'] for p in patients]
    lpid = len(pids)
    K = 10
    S = math.ceil(lpid/K)
    for i in range(0, lpid, S):
        test_p = patients[i:i+S]
        train_p = patients[:i] + patients[i+S:]

        linear_svc = SVC(kernel='linear')
        linear_svc.fit([p['P'] for p in train_p], [p['sick'] for p in train_p])
        is_hit = linear_svc.predict([p['P'] for p in test_p]) == [p['sick'] for p in test_p]
        # print('%3d: %s' % (i, is_hit))
        hit += np.sum(is_hit)
    print('Accuracy K-fold with K=%d accuracy=%.2lf' % (K, hit/len(pids)))
