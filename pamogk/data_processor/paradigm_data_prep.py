#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import collections
import pdb

from .. import config
from ..data_processor import rnaseq_processor as rp, synapse_rppa_processor as rpp
from ..gene_mapper import uniprot_mapper
from ..lib.sutils import *

parser = argparse.ArgumentParser(description='Run PAMOGK-mut algorithms on pathways')
parser.add_argument('--rs-patient-data', '-rs', metavar='file-path', dest='rnaseq_patient_data', type=Path,
                    help='rnaseq pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--rp-patient-data', '-rp', metavar='file-path', dest='rppa_patient_data', type=Path,
                    help='rppa pathway ID list', default=config.DATA_DIR / 'kirc_data/kirc_rppa_data')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=Path,
                    help='som mut pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/kirc_somatic_mutation_data.csv')
args = parser.parse_args()
print_args(args)


class ParadigmDataPrep(object):
    def __init__(self):
        self.exp_data_dir = config.DATA_DIR / 'paradigm'

        safe_create_dir(self.exp_data_dir)
        # change log and create log file
        change_log_path(self.exp_data_dir / 'run.log')
        log('exp_data_dir:', self.exp_data_dir)

    @timeit
    def read_rnaseq_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp, gene_name_map = rp.process(args.rnaseq_patient_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids
        ent_ids = gene_exp.index.values  # gene entrez ids
        return gene_exp.values, pat_ids, ent_ids

    @timeit
    def read_rppa_data(self):
        # Real Data #
        # process RNA-seq expression data

        gene_exp = rpp.process(args.rppa_patient_data)

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
        intersect_loc = self.exp_data_dir / 'patients.csv'
        with open(intersect_loc, 'w') as f:
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
    def save_data(self, fname, data, pat_ids, uni_ids):
        with (self.exp_data_dir / fname).open('w') as f:
            f.write('\t'.join(['gene_id', *[uids[0] for uids in uni_ids]]))
            for i, pid in enumerate(pat_ids):
                f.write('\n' + '\t'.join([str(v) for v in ['-'.join(pid.split('-')[:3]), *data[i][:]]]))


def main():
    exp = ParadigmDataPrep()

    # Patient part
    # RnaSeq Data
    rs_GE, rs_pat_ids, rs_ent_ids = exp.read_rnaseq_data()

    # Rppa Data
    rp_GE, rp_pat_ids, rp_ent_ids = exp.read_rppa_data()

    # Somatic mutation data
    som_patients = exp.read_som_data()

    # Find intersect
    rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients = exp \
        .find_intersection_patients(rs_GE, rs_pat_ids, rp_GE, rp_pat_ids, som_patients)

    # Kernel part
    # RnaSeq Data
    rs_GE, rs_uni_ids = exp.preprocess_seq_patient_data(rs_GE, rs_ent_ids)
    exp.save_data('rnaseq.tsv', rs_GE.T, rs_pat_ids, rs_uni_ids)

    # Rppa Data
    rp_GE, rp_uni_ids = exp.preprocess_seq_patient_data(rp_GE, rp_ent_ids)

    # Somatic mutation data
    som_patients = exp.preprocess_som_patient_data(som_patients)

    pdb.set_trace()


if __name__ == '__main__':
    main()
