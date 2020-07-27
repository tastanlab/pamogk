#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

from pamogk import config
from pamogk.data_processor import rnaseq_processor as rp
from pamogk.gene_mapper import uniprot_mapper
from pamogk.lib.sutils import *

parser = argparse.ArgumentParser(description='Run rMKL-LLP')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=Path,
                    help='Patient data file (if relative searched under data folder)',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
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


def main(*nargs):
    global args
    if __name__ == '__main__': # if running directly use command line arguments
        args = parser.parse_args()
    else: # otherwise use user given arguments
        args = parser.parse_args(nargs)

    exp = Experiment1()

    GE, pat_ids, ent_ids = exp.read_data()

    GE, uni_ids = exp.preprocess_patient_data(GE, ent_ids)

    with open('../ids.txt', 'w') as f:
        f.write('\n'.join(pat_ids))

    import scipy.io

    scipy.io.savemat('../file.mat', {'data': GE})

    log('Finished')

if __name__ == '__main__':
    main()
