#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import copy
import pdb
from operator import itemgetter
import operator
import matplotlib.pyplot as plt
import pandas as pd
from kernels.lmkkmeans_train import lmkkmeans_train
from data_processor import node2vec_processor
from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper
from kernels import center_product_kernel
import time
import json
import os
import csv
import config
from lib.sutils import *
import argparse
from sklearn.cluster import KMeans

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=str, help='pathway ID list', default='data/kirc_data/kirc_somatic_mutation_data.csv')
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable Debug Mode')
parser.add_argument('--disable-cache', '-c', dest='cache', action='store_false', help='disables intermediate caches')
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float, help='Node2Vec feature space size', default=128)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed', default=False)
parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments', default=1000)
parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float, help='Surviving patient percentage in range [0, 1]', default=0.9)
parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]', default=0.4)

args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self):
        self.exp_data_dir = os.path.join(config.data_dir, 'stmk', self.__class__.__name__)
        safe_create_dir(self.exp_data_dir)
        change_log_path(os.path.join(self.exp_data_dir, 'log-run={}.log'.format(args.rid)))
        log('exp_data_dir:', self.exp_data_dir)

        self.exp_result_dir = os.path.join(config.root_dir, '..', 'results')
        safe_create_dir(self.exp_result_dir)

    @timeit
    def read_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        patients = {}
        with open(args.patient_data) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                ent_id = row['Entrez Gene ID']
                if pat_id not in patients: patients[pat_id] = set([ent_id])
                else: patients[pat_id].add(ent_id)

        return patients

    @timeit
    def preprocess_patient_data(self, patients):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        res = []
        num_empty = 0
        for pat_id, ent_ids in patients.items():
            # uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            # if there are any matches map them
            if len(uni_ids) > 0: res.append({
                'pat_id': pat_id,
                'mutated_nodes': uni_ids,
            })
            else: num_empty += 1

        log('removed patients:', num_empty)

        return res

    @timeit
    def read_pathways(self):
        # get all pathways
        return cx_pw.read_pathways()

    @timeit
    def get_node2vecs(self, all_pw_map):
        fpath = 'p={p}-q={q}-size={n2v_size}-is_directed={is_directed}.json'.format(**vars(args))
        fpath = os.path.join(self.exp_data_dir, fpath)
        res = {}
        # if exists load from disk
        if args.cache and os.path.exists(fpath):
            with open(fpath) as f:
                return json.load(f)
        # otherwise calculate
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            log('Calculating node2vec for {:3}/{} pw_id={}'.format(ind + 1, num_pw, pw_id), end='\r')
            res[pw_id] = node2vec_processor.process(pw_id, pw, args, lambda x: x.tolist())
        log()
        # store gene vectors
        with open(fpath, 'w') as f: json.dump(res, f)
        return res

    @timeit
    def process_gene_vec_map(self, gene_vec_map, all_pw_map):
        uni_to_vec = {}
        for pw_id, pw_genes in gene_vec_map.items():
            pw = all_pw_map[pw_id]
            for n, gene_vec in pw_genes.items():
                for gene_id in pw.nodes[int(n)]['uniprot-ids']:
                    uni_to_vec[gene_id] = gene_vec
        return uni_to_vec

    @timeit
    def create_kernels(self, patients, gene_vec_map, uni_to_vec):
        center_product_kernel.calculate_S_and_P1(patients, gene_vec_map, uni_to_vec)
        return center_product_kernel.CP_kernels(patients)

    @timeit
    def cluster(self, kernels):
        return lmkkmeans_train(kernels)

    def save_results(self, **kwargs):
        save_np_data(os.path.join(self.exp_result_dir, 'stmk-exp-1-run={}'.format(args.rid)), **kwargs)



exp = Experiment1()

patient_map = exp.read_data()

patients = exp.preprocess_patient_data(patient_map)

all_pw_map = exp.read_pathways()

gene_vec_map = exp.get_node2vecs(all_pw_map)

uni_to_vec = exp.process_gene_vec_map(gene_vec_map, all_pw_map)

k1, k2 = exp.create_kernels(patients, gene_vec_map, uni_to_vec)

labels, h_normalized = exp.cluster(np.stack((k1,k2)))

exp.save_results(labels=labels, h_normalized=h_normalized)

log('Finished')

