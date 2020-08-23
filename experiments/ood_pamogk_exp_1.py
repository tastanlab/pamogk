#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

import networkx as nx

import pamogk
import pamogk.kernels.pamogk
from pamogk import config
from pamogk import label_mapper
from pamogk.data_processor import rnaseq_processor as rp
from pamogk.gene_mapper import uniprot_mapper
from pamogk.kernels.lmkkmeans_train import lmkkmeans_train
from pamogk.lib.sutils import *
from pamogk.pathway_reader import cx_pathway_reader as cx_pw

parser = argparse.ArgumentParser(description='Run PAMOGK algorithms on pathways')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=Path,
                    help='Patient data file (if relative searched under data folder)',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--disable-cache', '-c', dest='cache', action='store_false', help='disables intermediate caches')
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)

args = {}


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

        param_suffix = f'-label={label}-smoothing_alpha={smoothing_alpha}-norm={normalization}'
        exp_subdir = self.__class__.__name__ + param_suffix

        self.exp_data_dir = config.DATA_DIR / 'pamogk' / exp_subdir
        safe_create_dir(self.exp_data_dir)

        self.exp_result_dir = config.ROOT_DIR.parent / 'results'
        safe_create_dir(self.exp_result_dir)
        # change log and create log file
        change_log_path(self.exp_data_dir / f'log-run={args.rid}.log')
        log('exp_data_dir:', self.exp_data_dir)

        self.get_pw_path = lambda pw_id: self.exp_data_dir / f'pamogk-over-under-expressed-pw_id={pw_id}.gpickle'

    @timeit
    def read_data(self):
        # Real Data #
        # process RNA-seq expression data
        gene_exp, gene_name_map = rp.process(args.patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        log('num_pat:', pat_ids.shape)
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
            return self.get_pw_path(pw_id).exists()

        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    @timeit
    def restore_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_pw_path(pw_id)
            logr(f'Loading over/under expressed data {ind + 1:3}/{num_pw} pw={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_pw_path(pw_id)
            logr(f'Saving over/under expressed data {ind + 1:3}/{num_pw} pw={pw_id}')
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
            log(f'Checking patient for over-expressed  {ind + 1:4}/{num_pat} pid={pid}')
            gene_ind = (ge[..., pat_ids == pid] == 1).flatten()  # over expressed genes
            genes = uni_ids[gene_ind]  # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

            log(f'Checking patient for under-expressed {ind + 1:4}/{num_pat} pid={pid}')
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
        kms_path = self.exp_data_dir / 'kms.npz'
        if args.cache and kms_path.exists():
            return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            over_exp_kms[ind] = pamogk.kernels.pamogk.kernel(pat_ids, pw, label_key='label-oe', alpha=self.smoothing_alpha,
                                                             normalization=self.normalization)
            logr(f'Calculating oe pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            under_exp_kms[ind] = pamogk.kernels.pamogk.kernel(pat_ids, pw, label_key='label-ue', alpha=self.smoothing_alpha,
                                                              normalization=self.normalization)
            logr(f'Calculating ue pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster(self, kernels):
        return lmkkmeans_train(kernels)

    def save_results(self, **kwargs):
        np_save_npz(self.exp_result_dir / f'pamogk-exp-1-run={args.rid}', **kwargs)


def main(*nargs):
    global args
    if __name__ == '__main__':  # if running directly use command line arguments
        args = parser.parse_args()
    else:  # otherwise use user given arguments
        args = parser.parse_args(nargs)

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
