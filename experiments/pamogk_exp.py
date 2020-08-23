#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import collections

import mkkm_mr
import networkx as nx
from sklearn.cluster import KMeans, SpectralClustering
from snf_simple import SNF

from pamogk import config
from pamogk import label_mapper
from pamogk.data_processor import rnaseq_processor as rp, synapse_rppa_processor as rpp
from pamogk.gene_mapper import uniprot_mapper
from pamogk.kernels.lmkkmeans_train import lmkkmeans_train
from pamogk.kernels.pamogk import kernel
from pamogk.lib.sutils import *
from pamogk.pathway_reader import cx_pathway_reader as cx_pw
# see https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
from pamogk.result_processor.label_analysis import LabelAnalysis

# import sys
# sys.path.insert(0, '/Users/fma/dev/bilkent/research/snf')
# sys.path.insert(0, '/Users/fma/dev/bilkent/research/mkkm-mr')

parser = argparse.ArgumentParser(description='Run PAMOGK-mut algorithms on pathways')
parser.add_argument('--run-id', '-rid', metavar='run-id', dest='run_id', type=str, help='Unique Run ID')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=str2path,
                    help='rnaseq pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=str2path,
                    help='rppa pathway ID list', default=config.DATA_DIR / 'kirc_data/kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=str2path,
                    help='som mut pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/kirc_somatic_mutation_data.csv')
parser.add_argument('--label', '-m', metavar='label', dest='label', type=str, default='th196',
                    help='Label value that will be smoothed')
# used values: [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
parser.add_argument('--smoothing-alpha', '-a', metavar='alpha', dest='smoothing_alpha', type=float, default=0.01,
                    help='Smoothing alpha in range of 0-1')
parser.add_argument('--drop-percent', '-p', metavar='drop-percent', dest='drop_percent', type=int, default=1,
                    help='Drop percentage in range of 0-100')
parser.add_argument('--threshold', '-t', metavar='threshold', dest='threshold', type=float, default=1.96,
                    help='Cut off threshold')
parser.add_argument('--continuous', '-c', metavar='bool', dest='continuous', type=str2bool, default=True,
                    help='Whether to produce continuous values for under/over expressed')
parser.add_argument('--normalize-kernels', '-nk', dest='kernel_normalization', type=str2bool, default=True,
                    help='Kernel Normalization')
args = {}


class Experiment1(object):
    def __init__(self, args):
        """
        Parameters
        ----------
        args:
            arguments
        """
        self.args = args
        self.label = args.label
        self.smoothing_alpha = args.smoothing_alpha
        self.kernel_normalization = args.kernel_normalization
        self.drop_percent = args.drop_percent
        self.threshold = args.threshold
        self.log2_lambdas = list(range(-15, 16, 3))

        # these are kernel related params

        # each experiment may have different methods to build kernels
        exp_subdir = f'{Path(__file__).stem}-{self.__class__.__name__}'
        param_dir = f'label={self.label}-smoothing_alpha={self.smoothing_alpha}-kr_norm={self.kernel_normalization}'
        run_suffix = ''
        if self.args.run_id is not None:
            run_suffix = f'-run={self.args.run_id}'

        self.data_dir = config.DATA_DIR / 'pamogk_kirc' / exp_subdir / param_dir
        self.result_dir = self.data_dir / ('results' + run_suffix)
        self.kernel_dir = self.data_dir / 'kernels'

        self.label_analyzer = None

        # this will create with all roots
        safe_create_dir(self.result_dir)
        safe_create_dir(self.kernel_dir)
        # change log and create log file
        change_log_path(self.data_dir / 'run.log')
        log('exp_data_dir:', self.data_dir)

        self.get_rnaseq_pw_path = lambda \
                pw_id: self.kernel_dir / f'rnaseq-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_rppa_pw_path = lambda \
                pw_id: self.kernel_dir / f'rppa-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_som_pw_path = lambda \
                pw_id: self.kernel_dir / f'pamogk-som-expressed-pw_id={pw_id}.gpickle'

    @timeit
    def read_rnaseq_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp, gene_name_map = rp.process(self.args.rnaseq_patient_data, self.args.continuous, self.args.threshold)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp = rpp.process(self.args.rppa_patient_data, self.args.continuous, self.args.threshold)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_som_data(self):
        """
        Returns
        -------
        mapping of patient to mutations by entrez ids
        """
        # Real Data #
        # process RNA-seq expression data
        patients = {}
        with open(config.get_safe_data_file(self.args.som_patient_data)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                ent_id = row['Entrez Gene ID']
                if pat_id not in patients:
                    patients[pat_id] = {ent_id}
                else:
                    patients[pat_id].add(ent_id)

        return collections.OrderedDict(sorted(patients.items()))

    @timeit
    def find_intersection_patients(self, rs_GE, rs_pat, rp_GE, rp_pat, som_pat):
        rs_pat_list = simplify_pat_ids(rs_pat)
        rp_pat_list = simplify_pat_ids(rp_pat)
        som_pat_list = simplify_pat_ids(som_pat.keys())

        intersection_list = list(set(rs_pat_list).intersection(rp_pat_list, som_pat_list))
        intersection_list.sort()
        intersect_loc = self.data_dir / 'patients.csv'
        save_csv(intersect_loc, [[pid] for pid in intersection_list])

        def clean_patient_list_and_ge_data(patients, ge, whitelist):
            pat_list = simplify_pat_ids(patients)
            to_del = [idx for idx, value in enumerate(pat_list) if value not in whitelist]
            return np.delete(patients, to_del), np.delete(ge, to_del, axis=1)

        rs_pat, rs_GE = clean_patient_list_and_ge_data(rs_pat, rs_GE, intersection_list)
        rp_pat, rp_GE = clean_patient_list_and_ge_data(rp_pat, rp_GE, intersection_list)

        som_pat_deleted_list = [pid for pid in som_pat.keys() if pid not in intersection_list]

        for item in som_pat_deleted_list:
            som_pat.pop(item, None)

        return rs_GE, rs_pat, rp_GE, rp_pat, som_pat

    @timeit
    def preprocess_seq_patient_data(self, GE, all_ent_ids):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        found_ent_ids = [eid in ent2uni for eid in all_ent_ids]
        ent_ids = np.array([eid for eid in all_ent_ids if eid in ent2uni])
        uni_ids = np.array([ent2uni[eid] for eid in ent_ids], dtype=object)

        log('uni_ids:', len(uni_ids))
        log('miss_ent_ids:', len(all_ent_ids) - sum(found_ent_ids))

        # prune genes whose uniprot id is not found
        GE = GE[found_ent_ids]
        return GE, uni_ids

    @timeit
    def preprocess_som_patient_data(self, patients):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        res = []
        num_empty = 0
        for pat_id, ent_ids in patients.items():
            # uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            # if there are any matches map them
            res.append({
                'pat_id': pat_id,
                'mutated_nodes': uni_ids,
            })
        log('removed patients:', num_empty)

        return res

    @timeit
    def read_pathways(self):
        # get all pathways
        return cx_pw.read_pathways()

    def rnaseq_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_rnaseq_pw_path(pw_id).exists() for pw_id in all_pw_map])

    def rppa_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_rppa_pw_path(pw_id).exists() for pw_id in all_pw_map])

    def som_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_som_pw_path(pw_id).exists() for pw_id in all_pw_map])

    @timeit
    def restore_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr(f'Loading over/under rnaseq expressed data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rppa_pw_path(pw_id)
            logr(f'Loading over/under rppa expressed data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_som_pw_path(pw_id)
            logr(f'Loading somatic mutation data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr(f'Saving over/under rnaseq expressed data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rppa_pw_path(pw_id)
            logr(f'Saving over/under rppa expressed data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_som_pw_path(pw_id)
            logr(f'Saving somatic mutation data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def label_rnaseq_patient_genes(self, all_pw_map, pat_ids, GE, uni_ids):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        pat_ids: :obj:`list` of :obj:`str`
            list of patient ids
        GE: :obj:`numpy.ndarray`
            Gene expression data array in shape of genes by patients
        uni_ids: :obj:`numpy.ndarray`
            mapping from uniprot to gene
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rnaseq_pathways_save_valid(all_pw_map):
            return self.restore_rnaseq_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RNAseq Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            if self.args.continuous:
                gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
                logr(f'RNAseq Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
                label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
                label_mapper.mark_extra_label_on_pathways(f'oe-{self.label}', pid, all_pw_map, 'oe', self.threshold)

                logr(f'RNAseq Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
                label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
                label_mapper.mark_extra_label_on_pathways(f'ue-{self.label}', pid, all_pw_map, 'ue', self.threshold)
            else:
                logr(f'RNAseq Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
                gene_ind = (GE[..., pat_ids == pid] == 1).flatten()  # over expressed genes
                genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
                label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

                logr(f'RNAseq Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
                gene_ind = (GE[..., pat_ids == pid] == -1).flatten()  # under expressed genes
                genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
                label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)
        log()

        self.save_rnaseq_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def label_rppa_patient_genes(self, all_pw_map, pat_ids, GE, uni_ids):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        pat_ids: :obj:`list` of :obj:`str`
            list of patient ids
        GE: :obj:`numpy.ndarray`
            Gene expression data array in shape of genes by patients
        uni_ids: :obj:`numpy.ndarray`
            mapping from uniprot to gene
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rppa_pathways_save_valid(all_pw_map):
            return self.restore_rppa_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RPPA Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            if self.args.continuous:
                gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
                logr(f'RPPA Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
                label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
                label_mapper.mark_extra_label_on_pathways(f'oe-{self.label}', pid, all_pw_map, 'oe', self.threshold)

                logr(f'RPPA Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
                label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
                label_mapper.mark_extra_label_on_pathways(f'ue-{self.label}', pid, all_pw_map, 'ue', self.threshold)
            else:
                logr(f'RPPA Checking patient for rppa over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
                gene_ind = (GE[..., pat_ids == pid] == 1).flatten()  # over expressed genes
                genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
                label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

                logr(f'RPPA Checking patient for rppa under-expressed {ind + 1:4}/{num_pat} pid={pid}')
                gene_ind = (GE[..., pat_ids == pid] == -1).flatten()  # under expressed genes
                genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
                label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)
        log()

        self.save_rppa_pathways(all_pw_map)
        return all_pw_map

    def label_som_patient_genes(self, all_pw_map, patients):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        patients: :obj:`list`
            list of patients with mutation mappings
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.som_pathways_save_valid(all_pw_map):
            return self.restore_som_pathways(all_pw_map)

        num_pat = len(patients)
        # if there are missing ones calculate all of them
        log('Somatic mutation patient pathway labeling')
        for ind, patient in enumerate(patients):
            pid = patient['pat_id']
            genes = patient['mutated_nodes']  # get uniprot gene ids from indices
            genes = np.array([genes])
            logr(f'Checking patient for somatic mutation {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_label_on_pathways('som', pid, all_pw_map, genes, self.label)
        log()

        self.save_som_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_seq_kernels(self, all_pw_map, pat_ids, kms_file_name):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = self.kernel_dir / f'{kms_file_name}.npz'
        if kms_path.exists(): return np_load_data(kms_path, key='kms')
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            over_exp_kms[ind] = kernel(pat_ids, pw, label_key=f'label-oe-{self.label}', alpha=self.smoothing_alpha,
                                       normalization=self.kernel_normalization)
            logr(f'Calculating oe pathway kernel={kms_file_name} {ind + 1:4}/{num_pw} pw_id={pw_id}')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            under_exp_kms[ind] = kernel(pat_ids, pw, label_key=f'label-ue-{self.label}', alpha=self.smoothing_alpha,
                                        normalization=self.kernel_normalization)
            logr(f'Calculating ue pathway kernel={kms_file_name} {ind + 1:4}/{num_pw} pw_id={pw_id}')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def create_som_kernels(self, all_pw_map, patients):
        # experiment variables
        num_pat = len(patients)
        num_pw = len(all_pw_map)
        kms_path = self.kernel_dir / 'som-kms.npz'
        if kms_path.exists(): return np_load_data(kms_path, key='kms')
        # calculate kernel matrices for over expressed genes
        kms = np.zeros((num_pw, num_pat, num_pat))
        pat_ids = np.array([pat['pat_id'] for pat in patients])
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = kernel(pat_ids, pw, label_key='label-som', alpha=self.smoothing_alpha,
                              normalization=self.kernel_normalization)
            logr(f'Calculating som mut pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @staticmethod
    def kmeans_cluster(U, n_clusters):
        U_normalized = mkkm_mr.lib.normalize_unit_row(U)
        return KMeans(n_clusters=n_clusters, max_iter=100, n_init=50).fit_predict(U_normalized)

    def cluster_cont(self, kernels, n_clusters):
        snf_K = 20  # number of neighbors, usually (10~30)
        snf_t = 20  # number of iterations, usually (10~20)

        # SNF
        # W = snf_compute.snf(*kernels, K=snf_K, t=snf_t)
        W = SNF(kernels, K=snf_K, t=snf_t)

        # KMeans
        labels = self.kmeans_cluster(W, n_clusters)

        np_save_npz(self.result_dir / f'pamogk-snf-kmeans-k={n_clusters}', labels=labels)

        # Spectral
        labels = SpectralClustering(n_clusters, affinity='precomputed').fit_predict(W)
        np_save_npz(self.result_dir / f'pamogk-snf-spectral-k={n_clusters}', labels=labels)

        KH = mkkm_mr.lib.kernel_centralize(kernels)
        KH = mkkm_mr.lib.kernel_normalize(KH)
        num_ker = kernels.shape[0]
        gamma0 = np.ones((num_ker, 1)) / num_ker
        avgKer = mkkm_mr.lib.combine_kernels(KH, gamma0)

        H = mkkm_mr.lib.kernel_kmeans_iter(avgKer, n_clusters)
        labels = self.kmeans_cluster(H, n_clusters)
        np_save_npz(self.result_dir / f'pamogk-kmeans-k={n_clusters}.csv', labels=labels)

        # AAAI - 16 - MKKM-MR
        M = mkkm_mr.lib.calM(KH)
        lambdas = np.power(2., self.log2_lambdas)
        for log2_lambda, lambda_ in zip(self.log2_lambdas, lambdas):
            log(f'running for n_clusters={n_clusters} log2_lambda={log2_lambda}')
            [H, weights, obj] = mkkm_mr.mkkm_mr(KH, M, n_clusters, lambda_)
            labels = self.kmeans_cluster(H, n_clusters)
            out_file = self.result_dir / f'pamogk-mkkm-k={n_clusters}-log2_lambda={log2_lambda}'
            np_save_npz(out_file, labels=labels, weights=weights, obj=obj)

    def cluster_discrete(self, kernels, n_clusters):
        save_path = self.result_dir / f'labels_dropped={self.drop_percent}' / f'pamogk-all-lmkkmeans-k={n_clusters}'

        if save_path.exists():
            with np.load(save_path) as data:
                return data['labels', 'weights']

        labels, weights = lmkkmeans_train(kernels, cluster_count=n_clusters, iteration_count=5)
        ensure_file_dir(save_path)
        np_save_npz(f'{save_path}-weights', labels=labels, weights=weights)
        return labels, weights

    @timeit
    def cluster(self, kernels, n_clusters):
        if self.args.continuous:
            return self.cluster_cont(kernels, n_clusters)
        else:
            return self.cluster_discrete(kernels, n_clusters)

    @timeit
    def run(self):
        # Patient part
        # RnaSeq Data
        rs_GE, rs_pat_ids, rs_ent_ids = self.read_rnaseq_data()

        # Rppa Data
        rp_GE, rp_pat_ids, rp_ent_ids = self.read_rppa_data()

        # Somatic mutation data
        som_patients = self.read_som_data()

        # Find intersect
        rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients = self.find_intersection_patients(rs_GE, rs_pat_ids, rp_GE,
                                                                                             rp_pat_ids, som_patients)

        # Kernel part
        # RnaSeq Data
        rs_GE, rs_uni_ids = self.preprocess_seq_patient_data(rs_GE, rs_ent_ids)
        all_rs_pw_map = self.read_pathways()
        labeled_all_rs_pw_map = self.label_rnaseq_patient_genes(all_rs_pw_map, rs_pat_ids, rs_GE, rs_uni_ids)
        rs_kernels = self.create_seq_kernels(labeled_all_rs_pw_map, rs_pat_ids, 'rnaseq-kms')

        # Rppa Data
        rp_GE, rp_uni_ids = self.preprocess_seq_patient_data(rp_GE, rp_ent_ids)
        all_rp_pw_map = self.read_pathways()
        labeled_all_rp_pw_map = self.label_rppa_patient_genes(all_rp_pw_map, rp_pat_ids, rp_GE, rp_uni_ids)
        rp_kernels = self.create_seq_kernels(labeled_all_rp_pw_map, rp_pat_ids, 'rppa-kms')

        # Somatic mutation data
        som_patients = self.preprocess_som_patient_data(som_patients)
        all_som_pw_map = self.read_pathways()
        labeled_all_som_pw_map = self.label_som_patient_genes(all_som_pw_map, som_patients)
        som_kernels = self.create_som_kernels(labeled_all_som_pw_map, som_patients)

        kernels = np.concatenate((rs_kernels, rp_kernels, som_kernels))

        total = kernels.shape[1] * kernels.shape[2]
        limit = (self.drop_percent * total) / 100.0
        valid_kernels = kernels[np.count_nonzero(kernels, axis=(1, 2)) >= limit]

        log(f'kernel_count={kernels.shape[0]} valid_kernel_count={valid_kernels.shape[0]}')

        cluster_sizes = [2, 3, 4, 5]
        for k in cluster_sizes:
            log(f'Running clustering for k={k}')
            self.cluster(valid_kernels, k)

        self.label_analyzer = LabelAnalysis(results_dir=self.result_dir, methods=['mkkm', 'kmeans'],
                                            cluster_sizes=cluster_sizes, log2_lambdas=self.log2_lambdas)
        self.label_analyzer.run()


def create_experiment(*nargs):
    global args

    if __name__ == '__main__':  # if running directly use command line arguments
        args = parser.parse_args()
    else:  # otherwise use user given arguments
        args = parser.parse_args(nargs)

    print_args(args)

    return Experiment1(args)


if __name__ == '__main__':
    create_experiment().run()
