#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

import networkx as nx

import config
import label_mapper
import pamogk
from data_processor import rnaseq_processor as rp
from gene_mapper import uniprot_mapper
from kernels.lmkkmeans_train import lmkkmeans_train
from lib.sutils import *
from pathway_reader import cx_pathway_reader as cx_pw

parser = argparse.ArgumentParser(description='Run PAMOGK algorithms on pathways')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=str, help='pathway ID list',
                    default='../data/kirc_data/kirc_somatic_mutation_data.csv')
parser.add_argument('--disable-cache', '-c', dest='cache', action='store_false', help='disables intermediate caches')
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable Debug Mode')
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float,
                    help='Node2Vec feature space size', default=128)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed')
parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments',
                    default=1000)
parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float,
                    help='Surviving patient percentage in range [0, 1]', default=0.9)
parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]',
                    default=0.4)

args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self, label=1, smoothing_alpha=0, normalization=True):
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

        param_suffix = '-label={}-smoothing_alpha={}-norm={}'.format(label, smoothing_alpha, normalization)
        exp_subdir = self.__class__.__name__ + param_suffix

        self.exp_data_dir = os.path.join(config.data_dir, 'pamogk', exp_subdir)
        safe_create_dir(self.exp_data_dir)

        self.exp_result_dir = os.path.join(config.root_dir, '..', 'results')
        safe_create_dir(self.exp_result_dir)
        # change log and create log file
        change_log_path(os.path.join(self.exp_data_dir, 'log-run={}.log'.format(args.rid)))
        log('exp_data_dir:', self.exp_data_dir)

        data_file = 'pamogk-over-under-expressed'
        data_path = os.path.join(self.exp_data_dir, data_file)
        self.get_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(data_path, pw_id)

    @timeit
    def read_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        gene_exp, gene_name_map = rp.process(
            'data/kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        print('num_pat:', pat_ids.shape)
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def preprocess_patient_data(self, ge, all_ent_ids):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        found_ent_ids = [eid in ent2uni for eid in all_ent_ids]
        ent_ids = np.array([eid for eid in all_ent_ids if eid in ent2uni])
        uni_ids = np.array([ent2uni[eid] for eid in ent_ids])

        log('uni_ids:', len(uni_ids))
        log('miss_ent_ids:', len(all_ent_ids) - sum(found_ent_ids))

        # prune genes whose uniprot id is not found
        ge = ge[found_ent_ids]
        return ge, uni_ids

    @timeit
    def read_pathways(self):
        # get all pathways
        return cx_pw.read_pathways()

    def pathways_save_valid(self, all_pw_map):
        def pw_exists(pw_id):
            return os.path.exists(self.get_pw_path(pw_id))

        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    @timeit
    def restore_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_pw_path(pw_id)
            log('Loading over/under expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id), end='\r')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_pw_path(pw_id)
            log('Saving over/under expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id), end='\r')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def label_patient_genes(self, all_pw_map, pat_ids, ge, uni_ids):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        pat_ids: :obj:`list` of :obj:`str`
            list of patient ids
        ge: :obj:`numpy.ndarray`
            Gene expression data array in shape of genes by patients
        uni_ids: :obj:`numpy.ndarray`
            Uniprot gene id mapping
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.pathways_save_valid(all_pw_map):
            return self.restore_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            log('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (ge[..., pat_ids == pid] == 1).flatten()  # over expressed genes
            genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

            log('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (ge[..., pat_ids == pid] == -1).flatten()  # under expressed genes
            genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)

        self.save_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_kernels(self, all_pw_map, pat_ids):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = os.path.join(self.exp_data_dir, 'kms.npz')
        if args.cache and os.path.exists(kms_path):
            return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            over_exp_kms[ind] = pamogk.kernel(pat_ids, pw, label_key='label-oe', alpha=self.smoothing_alpha,
                                             normalization=self.normalization)
            log('Calculating oe pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id), end='\r')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            under_exp_kms[ind] = pamogk.kernel(pat_ids, pw, label_key='label-ue', alpha=self.smoothing_alpha,
                                              normalization=self.normalization)
            log('Calculating ue pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id), end='\r')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster(self, kernels):
        return lmkkmeans_train(kernels)

    def save_results(self, **kwargs):
        save_np_data(os.path.join(self.exp_result_dir, 'pamogk-exp-1-run={}'.format(args.rid)), **kwargs)


def main():
    exp = Experiment1()

    ge, pat_ids, ent_ids = exp.read_data()

    ge, uni_ids = exp.preprocess_patient_data(ge, ent_ids)

    all_pw_map = exp.read_pathways()

    labeled_all_pw_map = exp.label_patient_genes(all_pw_map, pat_ids, ge, uni_ids)

    kernels = exp.create_kernels(labeled_all_pw_map, pat_ids)

    labels, h_normalized = exp.cluster(kernels)

    exp.save_results(labels=labels, h_normalized=h_normalized)

    log('Finished')


if __name__ == '__main__':
    main()
