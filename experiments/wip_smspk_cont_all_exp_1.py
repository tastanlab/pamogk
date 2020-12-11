#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import collections
import collections as coll
import random
import sys

import matlab.engine

from pamogk.result_processor.label_analysis import LabelAnalysis

sys.path.append('..')
import networkx as nx
from pamogk import config
from pamogk.kernels.pamogk import kernel
from pamogk import label_mapper
from pamogk.data_processor import rnaseq_processor as rp, synapse_rppa_processor as rpp
from pamogk.gene_mapper import uniprot_mapper
from pamogk.lib.sutils import *
from pamogk.pathway_reader import cx_pathway_reader as cx_pw

# see https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html


parser = argparse.ArgumentParser(description='Run SMSPK-mut algorithms on pathways')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=str2path,
                    # help='rnaseq pathway ID list', default=Path('../data/kirc_data/rnaseqv2'))
                    help='rnaseq pathway ID list',
                    default=config.DATA_DIR / 'kirc_data' / 'unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=str2path,
                    # help='rppa pathway ID list', default=Path('../data/kirc_data/rppa'))
                    help='rppa pathway ID list', default=config.DATA_DIR / 'kirc_data' / 'kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=str2path,
                    help='som mut pathway ID list',
                    default=config.DATA_DIR / 'kirc_data' / 'kirc_somatic_mutation_data.csv')
parser.add_argument('--cancer-type', '-ct', metavar='cancer-type', dest='cancer_type', type=str,
                    help='som mut pathway ID list', default='kirc')
args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self, label=1, smoothing_alpha=0.05, normalization=True, bootstrap=None, bs_count=None):
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
        self.bootstrap = bootstrap
        self.bs_count = bs_count

        param_suffix = '-label={}-smoothing_alpha={}-norm={}'.format(label, smoothing_alpha, normalization)
        exp_subdir = self.__class__.__name__ + param_suffix

        param_suffix0 = '-label={}-smoothing_alpha={}-norm={}'.format(label, 0, normalization)
        exp_subdir0 = self.__class__.__name__ + param_suffix0

        self.orig_bs_pats_loc = None
        if self.bootstrap is not None:
            self.exp_data_dir = config.DATA_DIR / f'smspk_cont_{args.cancer_type}_all_paper_bs{self.bootstrap}_{bs_count}' / exp_subdir
            self.exp0_data_dir = config.DATA_DIR / f'smspk_cont_{args.cancer_type}_all_paper' / exp_subdir0  # _bs'+str(self.bootstrap)+'_'+str(bs_count), exp_subdir0)
            self.orig_bs_pats_loc = config.DATA_DIR / f'smspk_{args.cancer_type}_all_paper_bs{self.bootstrap}_{bs_count}' / exp_subdir, 'patients.csv'
        else:
            self.exp_data_dir = config.DATA_DIR / f'smspk_cont_{args.cancer_type}_all_paper' / exp_subdir
            self.exp0_data_dir = config.DATA_DIR / f'smspk_cont_{args.cancer_type}_all_paper' / exp_subdir0

        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(self.exp_data_dir / 'logs')
        log('exp_data_dir:', self.exp_data_dir)

        self.get_rnaseq_pw_path = lambda \
                pw_id: self.exp0_data_dir / f'smspk-rnaseq-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_rppa_pw_path = lambda \
                pw_id: self.exp0_data_dir / f'smspk-rppa-over-under-expressed-pw_id={pw_id}.gpickle'
        self.get_som_pw_path = lambda pw_id: self.exp0_data_dir / f'smspk-som-expressed-pw_id={pw_id}.gpickle'

    @timeit
    def read_rnaseq_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        gene_exp, gene_name_map = rp.process(args.rnaseq_patient_data, is_continuous=True)
        gene_exp = gene_exp.reindex(sorted(gene_exp.columns), axis=1)
        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        ### Real Data ###
        # process RNA-seq expression data

        gene_exp = rpp.process(args.rppa_patient_data, is_continuous=True)
        gene_exp = gene_exp.reindex(sorted(gene_exp.columns), axis=1)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
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
                if pat_id not in patients:
                    patients[pat_id] = {ent_id}
                else:
                    patients[pat_id].add(ent_id)
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

        if self.bootstrap is not None and self.orig_bs_pats_loc.exists():
            intersection_list = []
            with open(str(self.orig_bs_pats_loc), "r") as f:
                reader = csv.reader(f, delimiter=",")
                row1 = next(reader)
                for idx, pt_id in enumerate(row1):
                    intersection_list.append(pt_id)
        else:
            intersection_list = list(self.find_intersection_lists(rs_pat_list, rp_pat_list, som_pat_list))
            if self.bootstrap is not None:
                random.seed(self.bs_count)
                intersection_list = np.array(random.sample(intersection_list, k=self.bootstrap))
            intersection_list.sort()

        intersect_loc = self.exp_data_dir / "patients.csv"
        with open(intersect_loc, "w") as f:
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
        return np.all([self.get_rnaseq_pw_path(pw_id).exists() for pw_id in all_pw_map])

    def rppa_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_rppa_pw_path(pw_id).exists() for pw_id in all_pw_map])

    def som_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_som_pw_path(pw_id).exists() for pw_id in all_pw_map])

    @timeit
    def restore_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = coll.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr('Loading over/under rnaseq expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = coll.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_rppa_pw_path(pw_id)
            logr('Loading over/under rppa expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = coll.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_som_pw_path(pw_id)
            logr('Loading somatic mutation data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_rnaseq_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rnaseq_pw_path(pw_id)
            logr('Saving over/under rnaseq expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_rppa_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_rppa_pw_path(pw_id)
            logr('Saving over/under rppa expressed data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_som_pw_path(pw_id)
            logr('Saving somatic mutation data {:3}/{} pw={}'.format(ind + 1, num_pw, pw_id))
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
        uni_ids: list
            list of uniprot ids
        '''
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rnaseq_pathways_save_valid(all_pw_map):
            return self.restore_rnaseq_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RnaSeq Over and under expressed patient pathway labeling')
        # my_pw_map = {'68fce8e7-6193-11e5-8ac5-06603eb7f303': all_pw_map['68fce8e7-6193-11e5-8ac5-06603eb7f303']}
        for ind, pid in enumerate(pat_ids):
            gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
            logr('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways_fast('oe', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
            logr('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways_fast('ue', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
            # log('Checking patient for over-under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways('abs', pid, all_pw_map, uni_ids, gene_vals)
            # label_mapper.mark_extra_label_on_pathways('oe-th2', pid, all_pw_map, 'oe', 2)
            # label_mapper.mark_extra_label_on_pathways('ue-th2', pid, all_pw_map, 'ue', 2)
            label_mapper.mark_extra_label_on_pathways('oe-th196', pid, all_pw_map, 'oe', 1.96)
            label_mapper.mark_extra_label_on_pathways('ue-th196', pid, all_pw_map, 'ue', 1.96)
            # label_mapper.mark_extra_label_on_pathways('oe-th1', pid, all_pw_map, 'oe', 1.0)
            # label_mapper.mark_extra_label_on_pathways('ue-th1', pid, all_pw_map, 'ue', 1.0)
            # label_mapper.mark_extra_label_on_pathways('oe-th098', pid, all_pw_map, 'oe', 0.98)
            # label_mapper.mark_extra_label_on_pathways('ue-th098', pid, all_pw_map, 'ue', 0.98)
            # label_mapper.mark_extra_label_on_pathways('oe-th468', pid, all_pw_map, 'oe', 4.68)
            # label_mapper.mark_extra_label_on_pathways('ue-th468', pid, all_pw_map, 'ue', 4.68)
            # label_mapper.mark_extra_label_on_pathways('onek', pid, all_pw_map, 'xxx', 2)
        log()
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
        uni_ids: list
            list of uniprot ids
        '''
        # check if we already stored all over/under expression pathway data if so restore them
        if self.rppa_pathways_save_valid(all_pw_map):
            return self.restore_rppa_pathways(all_pw_map)

        num_pat = pat_ids.shape[0]
        # if there are missing ones calculate all of them
        log('RPPA Over and under expressed patient pathway labeling')
        for ind, pid in enumerate(pat_ids):
            gene_vals = (GE[..., pat_ids == pid]).flatten()  # over expressed genes
            log('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways_fast('oe', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_cont_label_on_pathways('oe', pid, all_pw_map, uni_ids, gene_vals)
            log('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways_fast('ue', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_cont_label_on_pathways('ue', pid, all_pw_map, uni_ids, gene_vals)
            label_mapper.mark_extra_label_on_pathways('oe-th196', pid, all_pw_map, 'oe', 1.96)
            label_mapper.mark_extra_label_on_pathways('ue-th196', pid, all_pw_map, 'ue', 1.96)
            # log('Checking patient for over-under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            # label_mapper.mark_cont_label_on_pathways('abs', pid, all_pw_map, uni_ids, gene_vals)
            # label_mapper.mark_extra_label_on_pathways('oe-th2', pid, all_pw_map, 'oe', 2)
            # label_mapper.mark_extra_label_on_pathways('ue-th2', pid, all_pw_map, 'ue', 2)
            # label_mapper.mark_extra_label_on_pathways('onek', pid, all_pw_map, 'xxx', 2)

        self.save_rppa_pathways(all_pw_map)
        return all_pw_map

    def label_som_patient_genes(self, all_pw_map, patients):
        '''Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        patients: :obj:`list` of :obj:`str`
            list of patient ids
        '''
        # check if we already stored all over/under expression pathway data if so restore them
        if self.som_pathways_save_valid(all_pw_map):
            return self.restore_som_pathways(all_pw_map)

        num_pat = len(patients)
        # if there are missing ones calculate all of them
        log('Somatic mutation patient pathway labeling')
        for ind, patient in enumerate(patients):
            pid = patient["pat_id"]
            genes = patient["mutated_nodes"]  # get uniprot gene ids from indices
            genes = np.array([genes])
            log('Checking patient for somatic mutation {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
            label_mapper.mark_label_on_pathways('som', pid, all_pw_map, genes, self.label)

        self.save_som_pathways(all_pw_map)
        return all_pw_map

    @timeit
    def create_seq_kernels(self, all_pw_map, pat_ids, kms_file_name, selection=''):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = self.exp_data_dir / f'{kms_file_name}-kms.npz'
        kms_abs_path = self.exp_data_dir / f'{kms_file_name}-abs-kms.npz'
        kms_onek_path = self.exp_data_dir / f'{kms_file_name}-onek-kms.npz'
        kms_th2_path = self.exp_data_dir / f'{kms_file_name}-th2-kms.npz'
        kms_th1_path = self.exp_data_dir / f'{kms_file_name}-th1-kms.npz'
        kms_th196_path = self.exp_data_dir / f'{kms_file_name}-th196-kms.npz'
        kms_th098_path = self.exp_data_dir / f'{kms_file_name}-th098-kms.npz'
        kms_th468_path = self.exp_data_dir / f'{kms_file_name}-thbonf-kms.npz'
        kms_th355_path = self.exp_data_dir / f'{kms_file_name}-thbonf-kms.npz'
        if selection == 'oe-ue' and kms_path.exists(): return np_load_data(kms_path, key='kms')
        if selection == 'abs' and kms_abs_path.exists(): return np.load(kms_abs_path)['kms']
        if selection == 'onek' and kms_onek_path.exists(): return np.load(kms_onek_path)['kms']
        if selection == 'th2' and kms_th2_path.exists(): return np.load(kms_th2_path)['kms']
        if selection == 'th1' and kms_th1_path.exists(): return np.load(kms_th1_path)['kms']
        if selection == 'th196' and kms_th196_path.exists(): return np.load(kms_th196_path)['kms']
        if selection == 'th098' and kms_th098_path.exists(): return np.load(kms_th098_path)['kms']
        if selection == 'th468' and kms_th468_path.exists(): return np.load(kms_th468_path)['kms']
        if selection == 'th355' and kms_th355_path.exists(): return np.load(kms_th355_path)['kms']
        kms_save_path = self.exp_data_dir / f'{kms_file_name}-{selection}-kms.npz'
        if selection == 'th355' or selection == 'th468':
            kms_save_path = self.exp_data_dir / f'{kms_file_name}-thbonf-kms.npz'

        if selection == 'oe-ue':
            # calculate kernel matrices for over expressed genes
            over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
            for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
                over_exp_kms[ind] = kernel(pat_ids, pw, label_key='label-oe', alpha=self.smoothing_alpha,
                                           normalization=self.normalization)
                logr('Calculating oe pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pw, pw_id))
            log()

            # calculate kernel matrices for under expressed genes
            under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
            for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
                under_exp_kms[ind] = kernel(pat_ids, pw, label_key='label-ue', alpha=self.smoothing_alpha,
                                            normalization=self.normalization)
                logr('Calculating ue pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pw, pw_id))
            log()
            kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
            np_save_npz(kms_path, kms=kms)  # save kernels
            return kms
        elif selection == 'abs' or selection == 'onek':  # calculate kernel matrices for absolute expressed genes
            abs_exp_kms = np.zeros((num_pw, num_pat, num_pat))
            for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
                abs_exp_kms[ind] = kernel(pat_ids, pw, label_key='label-' + selection, alpha=self.smoothing_alpha,
                                          normalization=self.normalization)
                logr('Calculating abs pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pw, pw_id))
            log()
            if selection == 'abs':
                np.savez_compressed(kms_abs_path, kms=abs_exp_kms)  # save kernels
            elif selection == 'onek':
                np.savez_compressed(kms_onek_path, kms=abs_exp_kms)  # save kernels
            return abs_exp_kms
        else:
            if len(selection.split('-')) > 1 and selection.split('-')[1] == 'th196' and kms_th196_path.exists():
                kms_ext = np.load(kms_th196_path)['kms']
                if selection.split('-')[0] == 'oe':
                    np.savez_compressed(kms_save_path, kms=kms_ext[:165])  # save kernels
                    return kms_ext[:165]
                if selection.split('-')[0] == 'ue':
                    np.savez_compressed(kms_save_path, kms=kms_ext[165:])  # save kernels
                    return kms_ext[165:]
                return 0
            # calculate kernel matrices for over expressed genes
            log('Running for cancer type: ', selection)
            over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
            for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
                over_exp_kms[ind] = kernel(pat_ids, pw, label_key='label-oe-' + selection,
                                           alpha=self.smoothing_alpha, normalization=self.normalization)
                logr('Calculating oe pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pw, pw_id))
            log()

            # calculate kernel matrices for under expressed genes
            under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
            for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
                under_exp_kms[ind] = kernel(pat_ids, pw, label_key='label-ue-' + selection,
                                            alpha=self.smoothing_alpha, normalization=self.normalization)
                logr('Calculating ue pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pw, pw_id))
            log()
            kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
            np.savez_compressed(kms_save_path, kms=kms)  # save kernels
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
        pat_ids = np.array([pat["pat_id"] for pat in patients])
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = kernel(pat_ids, pw, label_key='label-som', alpha=self.smoothing_alpha,
                              normalization=self.normalization)
            logr('Calculating som mut pathway kernel {:4}/{} pw_id={}'.format(ind + 1, num_pat, pw_id))
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster(self, kernels, cluster, drop_percent, typ=''):
        log('Clustering with smspk-cont')
        # return
        typ_c = ''
        if typ != '':
            typ_c = typ + '_'
        # Cluster using Mkkm-MR
        kmeans_save_path = self.exp_data_dir / f"labels_{typ_c}dropped{drop_percent}" / f"smspk-kmeans-{cluster}lab"

        if kmeans_save_path.exists():
            print('mkkm-mr already calculated')
        else:
            matlab_folder = config.DATA_DIR.parent / "pamogk_matlab"
            npy_matlab_folder1 = matlab_folder / "npy-matlab"
            snf_matlab_folder = matlab_folder / "SNFmatlab"
            npy_matlab_folder2 = npy_matlab_folder1 / "npy-matlab"
            eval_folder = matlab_folder / "ClusteringEvaluation"
            eng = matlab.engine.start_matlab()
            eng.addpath(str(npy_matlab_folder1))
            eng.addpath(str(npy_matlab_folder2))
            eng.addpath(str(matlab_folder))
            eng.addpath(str(eval_folder))
            eng.addpath(str(snf_matlab_folder))
            eng.addpath(str(self.exp_data_dir))
            # sending input to the function
            # eng.smspk_clustering_drop_fnc(str(self.exp_data_dir), cluster, drop_percent, typ)
            eng.smspk_clustering_fnc(str(self.exp_data_dir), cluster, drop_percent, typ)
            log('MKKM-MR and K-Means done.')

        # Cluster using lmkkmeans
        if typ == 'abs':
            save_path = self.exp_data_dir / f"labels_abs_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        elif typ == 'th2':
            save_path = self.exp_data_dir / f"labels_th2_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        elif typ == 'th1':
            save_path = self.exp_data_dir / f"labels_th1_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        elif typ == 'th196':
            save_path = self.exp_data_dir / f"labels_th196_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        elif typ == 'som':
            save_path = self.exp_data_dir / f"labels_som_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        elif typ == 'thbonf':
            save_path = self.exp_data_dir / f"labels_thbonf_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        else:
            save_path = self.exp_data_dir / f"labels_dropped{drop_percent}" / f"smspk-lmkkmeans-{cluster}lab"
        numsample = kernels.shape[1]
        '''
        if save_path.exists():
            print('lmmkmeans already calculated')
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

            results = lmkkmeans_train(kernels,cluster_count=cluster,iteration_count=50)
            ensure_file_dir(save_path)
            weights = np.mean(results[2], axis=0)
            weights = np.stack((stayed, weights))
            weights_loc = save_path+"weights"
            np.savetxt(weights_loc, weights.T, delimiter=",")
            obj_loc = save_path+"obj"
            np.savetxt(obj_loc,results[1] , delimiter=",")
            np.savetxt(save_path,results[0].labels_)
        '''
        return 1

    @timeit
    def callback(self):
        myList = []
        for i in range(330):
            name = "smspk-kernels-brca/" + str(i)
            myList.append(np.loadtxt(name))
        return np.array(myList)


def main():  #
    a_list = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    i_list = [1]
    # i_list = [14,15,16,17,18,19,20]
    bs = 300
    # a_list = [0.3]#, bootstrap=300, bs_count=i
    i_list = ['None']
    for i in i_list:
        for a in a_list:
            if i == 'None':
                bs = None
                i = None
            exp = Experiment1(smoothing_alpha=a, bootstrap=bs, bs_count=i)

            # Patient part
            # RnaSeq Data
            rs_GE, rs_pat_ids, rs_ent_ids = exp.read_rnaseq_data()

            # Rppa Data
            rp_GE, rp_pat_ids, rp_ent_ids = exp.read_rppa_data()

            # Somatic mutation data
            som_patients = exp.read_som_data()

            # Find intersect
            rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients = exp.find_intersection_patients(rs_GE, rs_pat_ids,
                                                                                                rp_GE, rp_pat_ids,
                                                                                                som_patients)

            # Kernel part
            # RnaSeq Data
            rs_GE, rs_uni_ids = exp.preprocess_seq_patient_data(rs_GE, rs_ent_ids)
            all_rs_pw_map = exp.read_pathways()
            labeled_all_rs_pw_map = exp.label_rnaseq_patient_genes(all_rs_pw_map, rs_pat_ids, rs_GE, rs_uni_ids)

            rs_kernels = exp.create_seq_kernels(labeled_all_rs_pw_map, rs_pat_ids, "rnaseq", 'th196')

            # Rppa Data
            rp_GE, rp_uni_ids = exp.preprocess_seq_patient_data(rp_GE, rp_ent_ids)
            all_rp_pw_map = exp.read_pathways()
            labeled_all_rp_pw_map = exp.label_rppa_patient_genes(all_rp_pw_map, rp_pat_ids, rp_GE, rp_uni_ids)
            rp_kernels = exp.create_seq_kernels(labeled_all_rp_pw_map, rp_pat_ids, "rppa", 'th196')

            # Somatic mutation data
            som_patients = exp.preprocess_som_patient_data(som_patients)
            all_som_pw_map = exp.read_pathways()
            labeled_all_som_pw_map = exp.label_som_patient_genes(all_som_pw_map, som_patients)
            som_kernels = exp.create_som_kernels(labeled_all_som_pw_map, som_patients)

            # print(rs_kernels.shape)
            # print(rp_kernels.shape)
            # print(som_kernels.shape)

            all_kernels = np.concatenate((rs_kernels, rp_kernels, som_kernels))
            for k in [2, 3, 4, 5]:
                log(f'Experimenting for k={k}')
                exp.cluster(all_kernels, k, 1, typ='th196')


if __name__ == '__main__':
    main()
