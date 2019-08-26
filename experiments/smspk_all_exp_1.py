#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys
sys.path.append('..')
import networkx as nx
import smspk
import label_mapper
from data_processor import rnaseq_processor as rp
from data_processor import synapse_rppa_processor as rpp
from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper
from kernels.lmkkmeans_train import lmkkmeans_train
import collections
import time
import config
from lib.sutils import *
import argparse
import csv


parser = argparse.ArgumentParser(description='Run SMSPK-mut algorithms on pathways')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=str, help='rnaseq pathway ID list', default='../data/kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=str, help='rppa pathway ID list', default='../data/kirc_data/kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=str, help='som mut pathway ID list', default='../data/kirc_data/kirc_somatic_mutation_data.csv')
args = parser.parse_args()
log('Running args:', args)

class Experiment1(object):
    def __init__(self, label = 1, smoothing_alpha = 0.05, normalization = True, drop_percent = 0):
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
        self.drop_percent = drop_percent

        param_suffix = '-label={}-smoothing_alpha={}-norm={}'.format(label, smoothing_alpha, normalization)
        exp_subdir = self.__class__.__name__ + param_suffix


        self.exp_data_dir = os.path.join(config.data_dir, 'smspk_kirc_all', exp_subdir)

        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(os.path.join(self.exp_data_dir, 'logs'))
        log('exp_data_dir:', self.exp_data_dir)

        rnaseq_data_file = 'smspk-rnaseq-over-under-expressed'
        rnaseq_data_path = os.path.join(self.exp_data_dir, rnaseq_data_file);
        self.get_rnaseq_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(rnaseq_data_path, pw_id)

        rppa_data_file = 'smspk-rppa-over-under-expressed'
        rppa_data_path = os.path.join(self.exp_data_dir, rppa_data_file);
        self.get_rppa_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(rppa_data_path, pw_id)

        som_data_file = 'smspk-som-expressed'
        som_data_path = os.path.join(self.exp_data_dir, som_data_file);
        self.get_som_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(som_data_path, pw_id)


    @timeit
    def read_rnaseq_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        gene_exp, gene_name_map = rp.process(args.rnaseq_patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values # patient TCGA ids
        ent_ids = gene_exp.index.values # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        gene_exp = rpp.process(args.rppa_patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values # patient TCGA ids
        ent_ids = gene_exp.index.values # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_som_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        patients = {}
        with open(args.som_patient_data) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                ent_id = row['Entrez Gene ID']
                if pat_id not in patients: patients[pat_id] = set([ent_id])
                else: patients[pat_id].add(ent_id)
        patients = collections.OrderedDict(sorted(patients.items()))

        return patients

    def find_intersection_lists(self, list1, list2, list3):
        intersection_list = set(list1).intersection(list2, list3)
        return intersection_list

    @timeit
    def find_intersection_patients(self, rs_GE, rs_pat, rp_GE, rp_pat, som_pat):
        rs_pat_list = []
        for pat in rs_pat:
            new_id = "-".join(pat.split("-")[0:3])
            rs_pat_list.append(new_id)

        rp_pat_list = []
        for pat in rp_pat:
            new_id = "-".join(pat.split("-")[0:3])
            rp_pat_list.append(new_id)

        som_pat_list = []
        for pat in som_pat.keys():
            som_pat_list.append(pat)

        intersection_list = list(self.find_intersection_lists(rs_pat_list, rp_pat_list, som_pat_list))
        intersection_list.sort()
        intersect_loc = os.path.join(self.exp_data_dir,"patients.csv")
        with open(intersect_loc,"w") as f:
            kirc_int = list(intersection_list)
            writer = csv.writer(f)
            writer.writerow(kirc_int)

        rs_pat_deleted_list = []
        for idx, value in enumerate(rs_pat_list):
            if value not in intersection_list:
                rs_pat_deleted_list.append(idx)

        rs_pat = np.delete(rs_pat, rs_pat_deleted_list)
        rs_GE = np.delete(rs_GE, rs_pat_deleted_list, axis=1)


        rp_pat_deleted_list = []
        for idx, value in enumerate(rp_pat_list):
            if value not in intersection_list:
                rp_pat_deleted_list.append(idx)

        rp_pat = np.delete(rp_pat, rp_pat_deleted_list)
        rp_GE = np.delete(rp_GE, rp_pat_deleted_list, axis=1)

        som_pat_deleted_list = []
        for pat_id in som_pat.keys():
            if pat_id not in intersection_list:
                som_pat_deleted_list.append(pat_id)

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
        pw_exists = lambda pw_id: os.path.exists(self.get_rnaseq_pw_path(pw_id))
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    def rppa_pathways_save_valid(self, all_pw_map):
        pw_exists = lambda pw_id: os.path.exists(self.get_rppa_pw_path(pw_id))
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    def som_pathways_save_valid(self, all_pw_map):
        pw_exists = lambda pw_id: os.path.exists(self.get_som_pw_path(pw_id))
        return np.all([pw_exists(pw_id) for pw_id in all_pw_map])

    @timeit
    def restore_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rnaseq_pw_path(pw_id)
            log('Loading over/under rnaseq expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rppa_pw_path(pw_id)
            log('Loading over/under rppa expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = {}
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_som_pw_path(pw_id)
            log('Loading somatic mutation data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rnaseq_pw_path(pw_id)
            log('Saving over/under rnaseq expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rppa_pw_path(pw_id)
            log('Saving over/under rppa expressed data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_som_pw_path(pw_id)
            log('Saving somatic mutation data {:3}/{} pw={}'.format(ind+1, num_pw, pw_id), end='\r')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def label_rnaseq_patient_genes(self, all_pw_map, pat_ids, GE, uni_ids):
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
        if self.rnaseq_pathways_save_valid(all_pw_map):
            return self.restore_rnaseq_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RnaSeq Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            log('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == 1).flatten() # over expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

            log('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == -1).flatten() # under expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)

        self.save_rnaseq_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def label_rppa_patient_genes(self, all_pw_map, pat_ids, GE, uni_ids):
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
        if self.rppa_pathways_save_valid(all_pw_map):
            return self.restore_rppa_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RPPA Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            log('Checking patient for rppa over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == 1).flatten() # over expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, self.label)

            log('Checking patient for rppa under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            gene_ind = (GE[..., pat_ids == pid] == -1).flatten() # under expressed genes
            genes = uni_ids[gene_ind] # get uniprot gene ids from indices
            label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, self.label)

        self.save_rppa_pathways(all_pw_map)
        return all_pw_map

    def label_som_patient_genes(self, all_pw_map, patients):
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
        if self.som_pathways_save_valid(all_pw_map):
            return self.restore_som_pathways(all_pw_map)

        num_pat = len(patients)
        # if there are missing ones calculate all of them
        log('Somatic mutation patient pathway labeling')
        for ind, patient in enumerate(patients):
            pid = patient["pat_id"]
            genes = patient["mutated_nodes"] # get uniprot gene ids from indices
            genes = np.array([genes])
            log('Checking patient for somatic mutation {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            label_mapper.mark_label_on_pathways('som', pid, all_pw_map, genes, self.label)

        self.save_som_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_seq_kernels(self, all_pw_map, pat_ids,kms_file_name):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = os.path.join(self.exp_data_dir, kms_file_name+'.npz')
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
    def create_som_kernels(self, all_pw_map, patients):
        # experiment variables
        num_pat = len(patients)
        num_pw = len(all_pw_map)
        kms_path = os.path.join(self.exp_data_dir, 'som-kms.npz')
        if os.path.exists(kms_path): return np.load(kms_path)['kms']
        # calculate kernel matrices for over expressed genes
        kms = np.zeros((num_pw, num_pat, num_pat))
        pat_ids = np.array([pat["pat_id"] for pat in patients])
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = smspk.kernel(pat_ids, pw, label_key='label-som', alpha=self.smoothing_alpha,
                                    normalization=self.normalization)
            log('Calculating som mut pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id), end='\r')
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster(self, kernels,cluster,drop_percent):
        save_path = os.path.join(self.exp_data_dir,"labels_dropped"+str(drop_percent),"smspk-all-lmkkmeans-"+str(cluster)+"lab")
        numsample = kernels.shape[1]
        if os.path.exists(save_path):
            return np.load(save_path)
        else:
            dropped = []
            stayed = []
            deletion = []
            total = numsample*numsample
            limit = (drop_percent*total)/100.0
            for i in range(len(kernels)):
                if np.count_nonzero(kernels[i]) < limit:
                    dropped.append(i+1)
                    deletion.append(i)
                else:
                    stayed.append(i+1)
            kernels = np.delete(kernels, deletion,axis=0)

            results = lmkkmeans_train(kernels[0:2],cluster_count=cluster,iteration_count=5)
            directory = os.path.dirname(save_path)
            safe_create_dir(directory)
            weights = np.mean(results[2], axis=0)
            weights = np.stack((stayed[0:2], weights))
            weights_loc = save_path+"weights"
            np.savetxt(weights_loc, weights.T, delimiter=",")
            np.save(save_path,results[0].labels_)
        return results[0].labels_


    @timeit
    def callback(self):
        myList = []
        for i in range(330):
            name = "smspk-kernels-brca/"+str(i)
            myList.append(np.loadtxt(name))
        return np.array(myList)


def main():
    for a in [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        exp = Experiment1(smoothing_alpha=a, drop_percent=1)

        # Patient part
        # RnaSeq Data
        rs_GE, rs_pat_ids, rs_ent_ids = exp.read_rnaseq_data()

        # Rppa Data
        rp_GE, rp_pat_ids, rp_ent_ids = exp.read_rppa_data()

        # Somatic mutation data
        som_patients = exp.read_som_data()

        #Find intersect
        rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients = exp.find_intersection_patients(rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients)


        # Kernel part
        # RnaSeq Data
        rs_GE, rs_uni_ids = exp.preprocess_seq_patient_data(rs_GE, rs_ent_ids)
        all_rs_pw_map = exp.read_pathways()
        labeled_all_rs_pw_map = exp.label_rnaseq_patient_genes(all_rs_pw_map, rs_pat_ids, rs_GE, rs_uni_ids)
        rs_kernels = exp.create_seq_kernels(labeled_all_rs_pw_map, rs_pat_ids, "rnaseq-kms")

        # Rppa Data
        rp_GE, rp_uni_ids = exp.preprocess_seq_patient_data(rp_GE, rp_ent_ids)
        all_rp_pw_map = exp.read_pathways()
        labeled_all_rp_pw_map = exp.label_rppa_patient_genes(all_rp_pw_map, rp_pat_ids, rp_GE, rp_uni_ids)
        rp_kernels = exp.create_seq_kernels(labeled_all_rp_pw_map, rp_pat_ids, "rppa-kms")

        # Somatic mutation data
        som_patients = exp.preprocess_som_patient_data(som_patients)
        all_som_pw_map = exp.read_pathways()
        labeled_all_som_pw_map = exp.label_som_patient_genes(all_som_pw_map, som_patients)
        som_kernels = exp.create_som_kernels(labeled_all_som_pw_map, som_patients)

        all_kernels = np.concatenate((rs_kernels, rp_kernels, som_kernels))

        for i in [2,3,4,5]:
            exp.cluster(all_kernels,i, exp.drop_percent)



if __name__ == '__main__':
    main()