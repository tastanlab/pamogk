#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

import config
from data_processor import rnaseq_processor as rp
from gene_mapper import uniprot_mapper
from lib.sutils import *

parser = argparse.ArgumentParser(description='Run rMKL-LLP')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=str,
                    help='patient data path', default='../data/kirc_data/kirc_somatic_mutation_data.csv')
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

        self.exp_result_dir = os.path.join(config.root_dir, '..', 'results')
        safe_create_dir(self.exp_result_dir)
        # change log and create log file
        change_log_path(os.path.join(self.exp_data_dir, 'log-run={}.log'.format(args.rid)))
        log('exp_data_dir:', self.exp_data_dir)

        data_file = 'smspk-over-under-expressed'
        data_path = os.path.join(self.exp_data_dir, data_file);
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


exp = Experiment1()

GE, pat_ids, ent_ids = exp.read_data()

GE, uni_ids = exp.preprocess_patient_data(GE, ent_ids)

with open('../ids.txt', 'w') as f:
    f.write('\n'.join(pat_ids))

import scipy.io

scipy.io.savemat('../file.mat', {'data': GE})

log('Finished')
