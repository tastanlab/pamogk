#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import copy
import pdb
from operator import itemgetter
import operator
import matplotlib.pyplot as plt
import smspk
import label_mapper
from data_processor import rnaseq_processor as rp
from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper
from kernels.lmkkmeans_train import lmkkmeans_train
import time
import json
import os
import config
from lib.sutils import *

class Experiment1(object):
    def __init__(self, label = 1, smoothing_alpha = 0, normalization = True):
        '''
        Parameters
        ----------
        label: {1} str
            label for over/under expressed
        smoothing_alpha: {0}
            smoothing parameter for smoothing out mutations
        '''
        self.label = label
        self.smoothing_alpha = smoothing_alpha
        self.normalization = normalization

        param_suffix = '-label={}-smoothing_alpha={}-norm={}'.format(label, smoothing_alpha, normalization)
        exp_subdir = self.__class__.__name__ + param_suffix

        self.exp_data_dir = os.path.join(config.data_dir, 'smspk', exp_subdir)
        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(os.path.join(self.exp_data_dir, 'logs'))
        log('exp_data_dir:', self.exp_data_dir)

        data_file = 'smspk-over-under-expressed'
        data_path = os.path.join(self.exp_data_dir, data_file);
        self.get_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(data_path, pw_id)

    @timeit
    def read_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        gene_exp, gene_name_map = rp.process('data/kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values # patient TCGA ids
        ent_ids = gene_exp.index.values # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def preprocess_patient_data(self, GE, all_ent_ids):
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
    def read_pathways(self):
        # get all pathways
        return cx_pw.read_pathways()

    def pathways_save_valid(self, all_pw_map):
        pw_exists = lambda pw_id: os.path.exists(self.get_pw_path(pw_id))
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    @timeit
    def restore_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_pw_path(pw_id)
            log('Loading over/under expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_pw_path(pw_id)
            log('Saving over/under expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def label_patient_genes(self, all_pw_map, pat_ids, GE):
        '''Labels all patients with matching level of expression

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
        '''
        # check if we already stored all over/under expression pathway data if so restore them
        if self.pathways_save_valid(all_pw_map):
            return self.restore_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            log('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == 1).flatten() # over expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

            log('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == -1).flatten() # under expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)

        self.save_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_kernels(self, all_pw_map, pat_ids):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = os.path.join(self.exp_data_dir, 'kms.npz')
        if os.path.exists(kms_path): return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()): # for each pathway
            over_exp_kms[ind] = smspk.kernel(pat_ids, pw, label_key='label-oe', alpha=self.smoothing_alpha, normalization=self.normalization)
            log('Calculating oe pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id), end='\r')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()): # for each pathway
            under_exp_kms[ind] = smspk.kernel(pat_ids, pw, label_key='label-ue', alpha=self.smoothing_alpha, normalization=self.normalization)
            log('Calculating ue pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id), end='\r')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms]) # stack all kernels
        np.savez_compressed(kms_path, kms=kms) # save kernels

        return kms

    @timeit
    def cluster(self, kernels):
        return lmkkmeans_train(kernels)

exp = Experiment1()

GE, pat_ids, ent_ids = exp.read_data()

GE, uni_ids = exp.preprocess_patient_data(GE, ent_ids)

all_pw_map = exp.read_pathways()

labeled_all_pw_map = exp.label_patient_genes(all_pw_map, pat_ids, GE)

kernels = exp.create_kernels(labeled_all_pw_map, pat_ids)

results = exp.cluster(kernels)

pdb.set_trace()
