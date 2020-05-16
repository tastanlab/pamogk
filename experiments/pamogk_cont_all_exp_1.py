#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import collections

import matlab.engine

import config
import label_mapper
import pamogk
from pamogk import rnaseq_processor as rp, synapse_rppa_processor as rpp
from pamogk import uniprot_mapper
from pamogk import lmkkmeans_train
from pamogk import *
from pamogk import cx_pathway_reader as cx_pw

parser = argparse.ArgumentParser(description='Run PAMOGK-mut algorithms on pathways')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=Path,
                    help='rnaseq pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=Path,
                    help='rppa pathway ID list', default=config.DATA_DIR / 'kirc_data/kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=Path,
                    help='som mut pathway ID list', default=config.DATA_DIR / 'kirc_data/kirc_somatic_mutation_data.csv')
args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self, label='th196', smoothing_alpha=0.05, normalization=True, drop_percent=0, thold=1.96):
        """
        Parameters
        ----------
        label: {1} str
            label for over/under expressed
        smoothing_alpha: {0}
            smoothing parameter for smoothing out mutations
        """
        self.label = label
        self.smoothing_alpha = smoothing_alpha
        self.normalization = normalization
        self.drop_percent = drop_percent
        self.thold = thold

        param_suffix = f'-label={label}-smoothing_alpha={smoothing_alpha}-norm={normalization}'
        exp_subdir = self.__class__.__name__ + param_suffix

        self.exp_data_dir = config.DATA_DIR / 'pamogk_cont_kirc_all' / exp_subdir

        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(self.exp_data_dir / 'run.log')
        log('exp_data_dir:', self.exp_data_dir)

        self.get_rnaseq_pw_path = lambda \
                pw_id: self.exp_data_dir / f'pamogk-rnaseq-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_rppa_pw_path = lambda \
                pw_id: self.exp_data_dir / f'pamogk-rppa-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_som_pw_path = lambda pw_id: self.exp_data_dir / f'pamogk-som-expressed-pw_id={pw_id}.gpickle'

    @timeit
    def read_rnaseq_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp, gene_name_map = rp.process_cont(args.rnaseq_patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp = rpp.process_cont(args.rppa_patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_som_data(self):
        # Real Data #
        # process RNA-seq expression data
        patients = {}
        with open(args.som_patient_data) as csvfile:
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
    def find_intersection_patients(self, rs_GE, rs_pat, rp_GE, rp_pat, som_pat):
        rs_pat_list = simplify_pat_ids(rs_pat)
        rp_pat_list = simplify_pat_ids(rp_pat)
        som_pat_list = simplify_pat_ids(som_pat.keys())

        intersection_list = list(set(rs_pat_list).intersection(rp_pat_list, som_pat_list))
        intersection_list.sort()
        intersect_loc = self.exp_data_dir / 'patients.csv'
        save_csv(intersect_loc, intersection_list)

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
        uni_ids = np.array([ent2uni[eid] for eid in ent_ids])

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
            '''
            if len(uni_ids) > 0: res.append({
                'pat_id': pat_id,
                'mutated_nodes': uni_ids,
            })
            else: num_empty += 1
            '''
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
        pw_exists = lambda pw_id: self.get_rnaseq_pw_path(pw_id).exists()
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    def rppa_pathways_save_valid(self, all_pw_map):
        pw_exists = lambda pw_id: self.get_rppa_pw_path(pw_id).exists()
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    def som_pathways_save_valid(self, all_pw_map):
        pw_exists = lambda pw_id: self.get_som_pw_path(pw_id).exists()
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    @timeit
    def restore_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr(f'Loading over/under rnaseq expressed data {ind + 1:3}/{num_pw} pw={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rppa_pw_path(pw_id)
            logr(f'Loading over/under rppa expressed data {ind + 1:3}/{num_pw} pw={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_som_pw_path(pw_id)
            logr(f'Loading somatic mutation data {ind + 1:3}/{num_pw} pw={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr(f'Saving over/under rnaseq expressed data {ind + 1:3}/{num_pw} pw={pw_id}')
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
            logr(f'Saving somatic mutation data {ind + 1:3}/{num_pw} pw={pw_id}')
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
        label: int, optional
            label that will be used for marking patients
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rnaseq_pathways_save_valid(all_pw_map):
            return self.restore_rnaseq_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RnaSeq Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
            log(f'Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_extra_label_on_pathways('oe-' + self.label, pid, all_pw_map, 'oe', thold=self.thold)

            log(f'Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_extra_label_on_pathways('ue-' + self.label, pid, all_pw_map, 'ue', thold=self.thold)

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
        label: int, optional
            label that will be used for marking patients
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rppa_pathways_save_valid(all_pw_map):
            return self.restore_rppa_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RPPA Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
            log(f'Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_extra_label_on_pathways('oe-' + self.label, pid, all_pw_map, 'oe', thold=self.thold)

            log(f'Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_extra_label_on_pathways('ue-' + self.label, pid, all_pw_map, 'ue', thold=self.thold)

        self.save_rppa_pathways(all_pw_map)
        return all_pw_map

    def label_som_patient_genes(self, all_pw_map, patients):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        patients: :obj:`list`
            list of patient
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
            log(f'Checking patient for somatic mutation {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_label_on_pathways('som', pid, all_pw_map, genes, self.label)

        self.save_som_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_seq_kernels(self, all_pw_map, pat_ids, kms_file_name):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = self.exp_data_dir / f'{kms_file_name}-{self.label}.npz'
        if kms_path.exists(): return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            over_exp_kms[ind] = pamogk.kernel(pat_ids, pw, label_key='label-oe-' + self.label,
                                              alpha=self.smoothing_alpha,
                                              normalization=self.normalization)
            logr(f'Calculating oe pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            under_exp_kms[ind] = pamogk.kernel(pat_ids, pw, label_key='label-ue-' + self.label,
                                               alpha=self.smoothing_alpha,
                                               normalization=self.normalization)
            logr(f'Calculating ue pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def create_som_kernels(self, all_pw_map, patients):
        # experiment variables
        num_pat = len(patients)
        num_pw = len(all_pw_map)
        kms_path = self.exp_data_dir / 'som-kms.npz'
        if kms_path.exists(): return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        kms = np.zeros((num_pw, num_pat, num_pat))
        pat_ids = np.array([pat['pat_id'] for pat in patients])
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = pamogk.kernel(pat_ids, pw, label_key='label-som', alpha=self.smoothing_alpha,
                                     normalization=self.normalization)
            logr(f'Calculating som mut pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster(self, kernels, cluster, drop_percent):
        log('Clustering with smspk-cont')
        # return
        typ_c = ''
        if self.label != '':
            typ_c = self.label + '_'
        # Cluster using Mkkm-MR
        mkkm_and_kkmeans_save_path = self.exp_data_dir / f'labels_{typ_c}dropped{str(drop_percent)}' / f'smspk-kmeans-{cluster}lab'

        if mkkm_and_kkmeans_save_path.exists():
            print('mkkm-mr and kk-means already calculated')
        else:
            matlab_folder = config.ROOT_DIR / 'pamogk_matlab'
            npy_matlab_folder1 = matlab_folder / 'npy-matlab'
            snf_matlab_folder = matlab_folder / 'SNFmatlab'
            npy_matlab_folder2 = npy_matlab_folder1 / 'npy-matlab'
            eval_folder = matlab_folder / 'ClusteringEvaluation'
            eng = matlab.engine.start_matlab()
            eng.addpath(npy_matlab_folder1)
            eng.addpath(npy_matlab_folder2)
            eng.addpath(matlab_folder)
            eng.addpath(eval_folder)
            eng.addpath(snf_matlab_folder)
            eng.addpath(self.exp_data_dir)
            eng.pamogk_clustering_fnc(self.exp_data_dir, cluster, drop_percent,
                                      self.label)  # sending input to the function
            log('MKKM-MR and K-Means done.')
        save_path = self.exp_data_dir / f'labels_{self.label}_dropped{drop_percent}' / f'pamogk-all-lmkkmeans-{cluster}lab'
        numsample = kernels.shape[1]
        if save_path.exists():
            return np.load(save_path)
        else:
            dropped = []
            stayed = []
            deletion = []
            total = numsample * numsample
            limit = (drop_percent * total) / 100.0
            for i in range(len(kernels)):
                if np.count_nonzero(kernels[i]) < limit:
                    dropped.append(i + 1)
                    deletion.append(i)
                else:
                    stayed.append(i + 1)
            kernels = np.delete(kernels, deletion, axis=0)

            results = lmkkmeans_train(kernels[0:2], cluster_count=cluster, iteration_count=5)
            ensure_file_dir(save_path)
            weights = np.mean(results[2], axis=0)
            weights = np.stack((stayed[0:2], weights))
            weights_loc = save_path + 'weights'
            np.savetxt(weights_loc, weights.T, delimiter=',')
            np.save(save_path, results[0].labels_)
        return results[0].labels_

    @timeit
    def callback(self):
        myList = []
        for i in range(330):
            name = f'pamogk-kernels-brca/{i}'
            myList.append(np.loadtxt(name))
        return np.array(myList)


def main():
    for a in [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        exp = Experiment1(label='th196', smoothing_alpha=a, drop_percent=1, thold=1.96)

        # Patient part
        # RnaSeq Data
        rs_GE, rs_pat_ids, rs_ent_ids = exp.read_rnaseq_data()

        # Rppa Data
        rp_GE, rp_pat_ids, rp_ent_ids = exp.read_rppa_data()

        # Somatic mutation data
        som_patients = exp.read_som_data()

        # Find intersect
        rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients = exp.find_intersection_patients(rs_GE, rs_pat_ids, rp_GE,
                                                                                            rp_pat_ids, som_patients)

        # Kernel part
        # RnaSeq Data
        rs_GE, rs_uni_ids = exp.preprocess_seq_patient_data(rs_GE, rs_ent_ids)
        all_rs_pw_map = exp.read_pathways()
        labeled_all_rs_pw_map = exp.label_rnaseq_patient_genes(all_rs_pw_map, rs_pat_ids, rs_GE, rs_uni_ids)
        rs_kernels = exp.create_seq_kernels(labeled_all_rs_pw_map, rs_pat_ids, 'rnaseq-kms')

        # Rppa Data
        rp_GE, rp_uni_ids = exp.preprocess_seq_patient_data(rp_GE, rp_ent_ids)
        all_rp_pw_map = exp.read_pathways()
        labeled_all_rp_pw_map = exp.label_rppa_patient_genes(all_rp_pw_map, rp_pat_ids, rp_GE, rp_uni_ids)
        rp_kernels = exp.create_seq_kernels(labeled_all_rp_pw_map, rp_pat_ids, 'rppa-kms')

        # Somatic mutation data
        som_patients = exp.preprocess_som_patient_data(som_patients)
        all_som_pw_map = exp.read_pathways()
        labeled_all_som_pw_map = exp.label_som_patient_genes(all_som_pw_map, som_patients)
        som_kernels = exp.create_som_kernels(labeled_all_som_pw_map, som_patients)

        all_kernels = np.concatenate((rs_kernels, rp_kernels, som_kernels))

        for i in [2, 3, 4, 5]:
            exp.cluster(all_kernels, i, exp.drop_percent)


if __name__ == '__main__':
    main()
