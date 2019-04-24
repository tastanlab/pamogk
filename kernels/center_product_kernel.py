import numpy as np
import pdb
from lib.sutils import *
import os,config,math
from sklearn.svm import SVC
from sklearn.metrics.pairwise import linear_kernel,rbf_kernel


def calculate_S_and_P(patients, gene_vectors):
    # calculate S (mutated gene vector set) and P (average mutataion point) vector
    for p in patients:
        genes = np.array([gene_vectors[n] for n in p['mutated_nodes']])
        p['S'] = genes
        p['P'] = np.average(genes, axis=0)

def calculate_S_and_P1(patients, gene_vectors, uni_to_vec):
    # calculate S (mutated gene vector set) and P (average mutataion point) vector
    for p in patients:
        genes = []
        for pw_genes in gene_vectors.values():
            for n in p['mutated_nodes']:
                genes.append(uni_to_vec[n])
        p['S'] = genes
        p['P'] = np.average(genes, axis=0)
    return patients

def CP_kernels(patients):
    vectors = np.array([p['P'] for p in patients])
    return linear_kernel(vectors), rbf_kernel(vectors)

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
    log('Accuracy Leave-One-Out accuracy=%.2lf' % (hit/len(pids)))

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
    log('Accuracy K-fold with K=%d accuracy=%.2lf' % (K, hit/len(pids)))
