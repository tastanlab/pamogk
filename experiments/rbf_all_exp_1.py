#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import collections
import csv
import math

from sklearn import preprocessing
from sklearn.metrics.pairwise import rbf_kernel

import config
from data_processor import rnaseq_processor as rp
from data_processor import synapse_rppa_processor as rpp
from kernels.lmkkmeans_train import lmkkmeans_train
from lib.sutils import *

parser = argparse.ArgumentParser(description='Run PAMOGK-mut algorithms on pathways')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=str,
                    help='RNAseq Patient data file (relative to data dir, otherwise use absolute path)',
                    default='kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=str,
                    help='RPPA Patient data file (relative to data dir, otherwise use absolute path)',
                    default='kirc_data/kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=str,
                    help='Somatic Patient data file (relative to data dir, otherwise use absolute path)',
                    default='kirc_data/kirc_somatic_mutation_data.csv')

args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self, label=1, gamma=None, normalization=True):
        """
        Parameters
        ----------
        label: {1} str
            label for over/under expressed
        gamma: {0}
            gamma parameter for rbf kernel
        """
        self.label = label
        self.gamma = gamma
        self.normalization = normalization

        param_suffix = '-label={}-gamma={}-norm={}'.format(label, gamma, normalization)
        exp_subdir = self.__class__.__name__ + param_suffix

        self.exp_data_dir = os.path.join(config.data_dir, 'rbf_kirc_all', exp_subdir)

        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(os.path.join(self.exp_data_dir, 'logs'))
        log('exp_data_dir:', self.exp_data_dir)

    @timeit
    def read_rnaseq_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        fpath = config.get_safe_data_file(args.rnaseq_patient_data)
        gene_exp, gene_name_map = rp.process(fpath)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        fpath = config.get_safe_data_file(args.rppa_patient_data)
        gene_exp = rpp.process(fpath)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_som_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        patients = {}
        fpath = config.get_safe_data_file(args.som_patient_data)
        with open(fpath) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                ent_id = row['Entrez Gene ID']
                if pat_id not in patients:
                    patients[pat_id] = {ent_id}
                else:
                    patients[pat_id].add(ent_id)
        patients = collections.OrderedDict(sorted(patients.items()))

        return patients

    @staticmethod
    def find_intersection_lists(list1, list2, list3):
        intersection_list = set(list1).intersection(list2, list3)
        return intersection_list

    @timeit
    def find_intersection_patients(self, rs_ge, rs_pat, rp_ge, rp_pat, som_pat):
        rs_pat_list = []
        for pat in rs_pat:
            new_id = '-'.join(pat.split('-')[0:3])
            rs_pat_list.append(new_id)

        rp_pat_list = []
        for pat in rp_pat:
            new_id = '-'.join(pat.split('-')[0:3])
            rp_pat_list.append(new_id)

        som_pat_list = []
        for pat in som_pat.keys():
            som_pat_list.append(pat)

        intersection_list = list(self.find_intersection_lists(rs_pat_list, rp_pat_list, som_pat_list))
        intersection_list.sort()
        intersect_loc = os.path.join(self.exp_data_dir, 'patients.csv')
        with open(intersect_loc, 'w') as f:
            kirc_int = list(intersection_list)
            writer = csv.writer(f)
            writer.writerow(kirc_int)

        rs_pat_deleted_list = []
        for idx, value in enumerate(rs_pat_list):
            if value not in intersection_list:
                rs_pat_deleted_list.append(idx)

        rs_pat = np.delete(rs_pat, rs_pat_deleted_list)
        rs_ge = np.delete(rs_ge, rs_pat_deleted_list, axis=1)

        rp_pat_deleted_list = []
        for idx, value in enumerate(rp_pat_list):
            if value not in intersection_list:
                rp_pat_deleted_list.append(idx)

        rp_pat = np.delete(rp_pat, rp_pat_deleted_list)
        rp_ge = np.delete(rp_ge, rp_pat_deleted_list, axis=1)

        som_pat_deleted_list = []
        for pat_id in som_pat.keys():
            if pat_id not in intersection_list:
                som_pat_deleted_list.append(pat_id)

        for item in som_pat_deleted_list:
            som_pat.pop(item, None)

        return rs_ge, rs_pat, rp_ge, rp_pat, som_pat

    @staticmethod
    def find_gammas(matrix):
        scaler = preprocessing.StandardScaler().fit(matrix)
        scaler.scale_ = np.std(matrix, axis=0, ddof=1)
        x = scaler.transform(matrix)
        m = x.shape[0]
        n = math.floor(0.5 * m)
        # index = np.floor(np.random.rand(n) * m)
        # index2 = np.floor(np.random.rand(n) * m)
        index = np.random.randint(0, high=m, size=n)
        index2 = np.random.randint(0, high=m, size=n)
        temp = x[index, :] - x[index2, :]
        dist = np.nansum(np.multiply(temp, temp), axis=1)
        dist = dist[np.nonzero(dist)]
        gammas = np.power(np.quantile(dist, [0.9, 0.5, 0.1]), -1)
        return gammas

    @staticmethod
    def sig_to_gamma(sigma):
        sigma2 = 2 * np.power(sigma, 2)
        return 1.0 / sigma2

    @staticmethod
    def find_sigma(data):
        size1 = data.shape[0]
        if size1 > 100:  # Choose 100 random samples from data if it contains more than 100 samples
            np.random.seed(34956)
            randoms = np.random.rand(size1)
            ind = np.argsort(randoms)
            data_med = data[ind[:100], :]
            size1 = 100
        else:
            data_med = data

        g = np.nansum(np.multiply(data_med, data_med), axis=1).reshape((size1, 1))
        q = np.tile(g, (1, size1))
        r = np.tile(g.T, (size1, 1))
        dists = q + r - np.matmul(data_med, data_med.T) * 2
        dists = dists - np.tril(dists)
        dists = dists.reshape((size1 ** 2, 1), order='F')
        dists = dists[dists > 0]
        sig = np.sqrt(0.5 * np.median(dists))
        return sig

    @timeit
    def create_seq_kernels(self, ge, kms_file_name):
        # experiment variables
        np.random.seed()
        save_over = os.path.join(self.exp_data_dir, kms_file_name + '-over')
        save_under = os.path.join(self.exp_data_dir, kms_file_name + '-under')
        save_over_under = os.path.join(self.exp_data_dir, kms_file_name + '-over-under')
        if os.path.exists(save_over + '.npy') and os.path.exists(save_under + '.npy') and \
                os.path.exists(save_over_under + '.npy'):
            kernel_over = np.load(save_over + '.npy')
            kernel_under = np.load(save_under + '.npy')
            kernel_over_under = np.load(save_over_under + '.npy')
            return kernel_over, kernel_under, kernel_over_under
        ge_t = ge.T
        ge_over = ge_t.copy()
        ge_under = ge_t.copy()
        ge_over_under = ge_t.copy()

        ge_over[ge_over == -1] = 0
        ge_under[ge_under == 1] = 0
        ge_under[ge_under == -1] = 1
        ge_over_under[ge_over_under == -1] = 1
        cur_gamma = self.gamma

        if self.gamma == 'median':
            gammas_over = self.find_gammas(ge_over)
            gammas_under = self.find_gammas(ge_under)
            gammas_over_under = self.find_gammas(ge_over_under)

            for med_gamma in gammas_over_under:
                kernel = rbf_kernel(ge_over_under, gamma=med_gamma)
                np.save('{}-gamma={:.2e}'.format(save_over_under, med_gamma), kernel)

            for med_gamma in gammas_over:
                kernel = rbf_kernel(ge_over_under, gamma=med_gamma)
                np.save('{}-gamma={:.2e}'.format(save_over, med_gamma), kernel)

            for med_gamma in gammas_under:
                kernel_under = rbf_kernel(ge_under, gamma=med_gamma)
                np.save('{}-gamma={:.2e}'.format(save_under, med_gamma), kernel_under)
            cur_gamma = None
        elif self.gamma == 'median2':
            med_gamma = self.sig_to_gamma(self.find_sigma(ge_over_under))
            kernel_over_under = rbf_kernel(ge_over_under, gamma=med_gamma)
            np.save('{}-gamma={:.2e}'.format(save_over_under, med_gamma), kernel_over_under)

            med_gamma = self.sig_to_gamma(self.find_sigma(ge_over))
            kernel_over = rbf_kernel(ge_over_under, gamma=med_gamma)
            np.save('{}-gamma={:.2e}'.format(save_over, med_gamma), kernel_over)

            med_gamma = self.sig_to_gamma(self.find_sigma(ge_under))
            kernel_under = rbf_kernel(ge_over_under, gamma=med_gamma)
            np.save('{}-gamma={:.2e}'.format(save_under, med_gamma), kernel_under)
            return kernel_over, kernel_under, kernel_over_under
        kernel_over = rbf_kernel(ge_over, gamma=cur_gamma)
        kernel_under = rbf_kernel(ge_under, gamma=cur_gamma)
        kernel_over_under = rbf_kernel(ge_over_under, gamma=cur_gamma)
        np.save(save_over, kernel_over)
        np.save(save_under, kernel_under)
        np.save(save_over_under, kernel_over_under)
        return kernel_over, kernel_under, kernel_over_under

    @timeit
    def create_som_kernels(self, patients):
        # experiment variables
        save_loc = os.path.join(self.exp_data_dir, 'som')
        if os.path.exists(save_loc + '.npy'):
            kernel = np.load(save_loc + '.npy')
            return kernel
        num_pat = len(patients)
        entrez_ids = []
        for patient, e_ids in patients.items():
            entrez_ids.extend(e_ids)
        e_ids = list(set(entrez_ids))
        num_gene = len(e_ids)
        feature_matrix = np.zeros((num_pat, num_gene))

        for col_id, e_id in enumerate(e_ids):
            for row_id, (patient, e_ids) in enumerate(patients.items()):
                if e_id in e_ids:
                    feature_matrix[row_id, col_id] = 1

        cur_gamma = self.gamma

        if self.gamma == 'median':

            gammas = self.find_gammas(feature_matrix)
            for med_gamma in gammas:
                kernel = rbf_kernel(feature_matrix, gamma=med_gamma)
                f_gamma = '{:.2e}'.format(med_gamma)
                np.save('{}-gamma={}'.format(save_loc, f_gamma), kernel)
            cur_gamma = None
        elif self.gamma == 'median2':
            sigma = self.find_sigma(feature_matrix)
            med_gamma = self.sig_to_gamma(sigma)
            kernel = rbf_kernel(feature_matrix, gamma=med_gamma)
            f_gamma = '{:.2e}'.format(med_gamma)
            np.save('{}-gamma={}'.format(save_loc, f_gamma), kernel)
            return kernel

        kernel = rbf_kernel(feature_matrix, gamma=cur_gamma)
        np.save(save_loc, kernel)

        return kernel

    @timeit
    def cluster(self, kernels, cluster, run_type):
        save_path = os.path.join(self.exp_data_dir, 'labels', 'rbf-all-{}-lmkkmeans-{}lab'.format(run_type, cluster))
        if os.path.exists(save_path):
            return np.load(save_path)
        else:
            labels, _ = lmkkmeans_train(kernels, cluster_count=cluster, iteration_count=5)
            directory = os.path.dirname(save_path)
            safe_create_dir(directory)
            np.save(save_path, labels)
        return labels

    @timeit
    def callback(self):
        return np.array([np.loadtxt('pamogk-kernels-brca/{}'.format(i)) for i in range(330)])


def main():
    gammas = ['median2']
    for g in gammas:
        exp = Experiment1(gamma=g)

        # Patient part
        # RnaSeq Data
        rs_ge, rs_pat_ids, rs_ent_ids = exp.read_rnaseq_data()

        # Rppa Data
        rp_ge, rp_pat_ids, rp_ent_ids = exp.read_rppa_data()

        # Somatic mutation data
        som_patients = exp.read_som_data()

        # Find intersect
        rs_ge, rs_pat_ids, rp_ge, rp_pat_ids, som_patients = exp \
            .find_intersection_patients(rs_ge, rs_pat_ids, rp_ge, rp_pat_ids, som_patients)

        # Kernel part
        # RnaSeq Data
        rs_kernels = exp.create_seq_kernels(rs_ge, 'rnaseq')

        # Rppa Data
        rp_kernels = exp.create_seq_kernels(rp_ge, 'rppa')

        # Somatic mutation data
        som_kernels = exp.create_som_kernels(som_patients)

        all_kernels1 = np.stack((rs_kernels[2], rp_kernels[2], som_kernels))
        # all_kernels2 = np.stack((rs_kernels[0],rs_kernels[1], rp_kernels[0],rp_kernels[1], som_kernels))

        for i in [2, 3, 4, 5]:
            exp.cluster(all_kernels1, i, '3ker')
            # exp.cluster(all_kernels2,i,'5ker')


if __name__ == '__main__':
    main()
