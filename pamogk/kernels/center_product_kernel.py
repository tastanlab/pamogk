import math

from sklearn.metrics.pairwise import linear_kernel
from sklearn.svm import SVC

from ..lib.sutils import *


def calculate_S_and_P(patients, gene_vectors, uni_to_vec):
    # calculate S (mutated gene vector set) and P (average mutataion point) vector
    kms_map = {}
    for p in patients:
        genes = []
        p['S'] = {}
        p['P'] = {}
        for pw_id, pw_genes in gene_vectors.items():
            for n in p['mutated_nodes']:
                genes.append(uni_to_vec[n])
            P = np.average(genes, axis=0)
            p['S'] = {pw_id: genes}
            p['P'] = {pw_id: P}
            if pw_id not in kms_map:
                kms_map[pw_id] = [P]
            else:
                kms_map[pw_id].append(P)
    # convert list to array
    for pw_id in kms_map:
        kms_map[pw_id] = np.vstack(kms_map[pw_id])
    return kms_map


def CP_kernels(kms, method=linear_kernel):
    return map(method, kms.values())


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
    log('Accuracy Leave-One-Out accuracy=%.2lf' % (hit / len(pids)))

    hit = 0
    pids = [p['pid'] for p in patients]
    lpid = len(pids)
    K = 10
    S = math.ceil(lpid / K)
    for i in range(0, lpid, S):
        test_p = patients[i:i + S]
        train_p = patients[:i] + patients[i + S:]

        linear_svc = SVC(kernel='linear')
        linear_svc.fit([p['P'] for p in train_p], [p['sick'] for p in train_p])
        is_hit = linear_svc.predict([p['P'] for p in test_p]) == [p['sick'] for p in test_p]
        # print('%3d: %s' % (i, is_hit))
        hit += np.sum(is_hit)
    log('Accuracy K-fold with K=%d accuracy=%.2lf' % (K, hit / len(pids)))
