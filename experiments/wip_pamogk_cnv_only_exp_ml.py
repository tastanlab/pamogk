#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import collections

import matlab.engine
import networkx as nx
import scipy.io

from pamogk import config
from pamogk import label_mapper
from pamogk.gene_mapper import uniprot_mapper
from pamogk.kernels.pamogk import kernel
from pamogk.lib.sutils import *
from pamogk.pathway_reader import cx_pathway_reader as cx_pw
# see https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
from pamogk.result_processor.label_analysis_new import LabelAnalysis

# import sys
# sys.path.insert(0, '/Users/fma/dev/bilkent/research/snf')
# sys.path.insert(0, '/Users/fma/dev/bilkent/research/mkkm-mr')

parser = argparse.ArgumentParser(description='Run PAMOGK-mut algorithms on pathways')
parser.add_argument('--run-id', '-rid', metavar='run-id', dest='run_id', type=str, help='Unique Run ID')
parser.add_argument('--som-patient-data', '-s', metavar='file-path', dest='som_patient_data', type=str2path,
                    help='som mut gene ID list',
                    default=config.DATA_DIR / 'kirc_data' / 'kirc_somatic_mutation_data.csv')
parser.add_argument('--cnv-patient-data', metavar='file-path', dest='cnv_patient_data', type=str2path,
                    help='CNV data with uniprot gene ids already mapped',
                    default=config.DATA_DIR / 'copy-number-variation' / 'KIRC.focal_score_by_genes-tcga_uniprot.csv')
parser.add_argument('--label', '-m', metavar='label', dest='label', type=str, default='th196',
                    help='Label value that will be smoothed')
# used values: [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
parser.add_argument('--smoothing-alpha', '-a', metavar='alpha', dest='smoothing_alpha', type=float, default=0.01,
                    help='Smoothing alpha in range of 0-1')
parser.add_argument('--drop-percent', '-p', metavar='drop-percent', dest='drop_percent', type=int, default=1,
                    help='Drop percentage in range of 0-100')
parser.add_argument('--threshold', '-t', metavar='threshold', dest='threshold', type=float, default=1.96,
                    help='Cut off threshold')
parser.add_argument('--continuous', '-c', metavar='bool', dest='continuous', type=str2bool, default=True,
                    help='Whether to produce continuous values for under/over expressed')
parser.add_argument('--normalize-kernels', '-nk', dest='kernel_normalization', type=str2bool, default=True,
                    help='Kernel Normalization')
args = {}


class Experiment1(object):
    def __init__(self, args):
        """
        Parameters
        ----------
        args:
            arguments
        """
        self.args = args
        self.label = args.label
        self.smoothing_alpha = args.smoothing_alpha
        self.kernel_normalization = args.kernel_normalization
        self.drop_percent = args.drop_percent
        self.threshold = args.threshold
        self.log2_lambdas = list(range(-15, 16, 3))

        # these are kernel related params

        # each experiment may have different methods to build kernels
        exp_subdir = f'{Path(__file__).stem}-{self.__class__.__name__}'
        param_dir = f'label={self.label}-smoothing_alpha={self.smoothing_alpha}-kr_norm={self.kernel_normalization}'
        run_suffix = ''
        if self.args.run_id is not None:
            run_suffix = f'-run={self.args.run_id}'

        self.data_dir = config.DATA_DIR / 'pamogk_kirc' / exp_subdir / param_dir
        self.result_dir = self.data_dir / ('results' + run_suffix)
        self.kernel_dir = self.data_dir / 'kernels'

        self.label_analyzer = None

        # this will create with all roots
        safe_create_dir(self.result_dir)
        safe_create_dir(self.kernel_dir)
        # change log and create log file
        change_log_path(self.data_dir / 'run.log')
        log('exp_data_dir:', self.data_dir)

        self.get_som_pw_path = lambda \
                pw_id: self.kernel_dir / f'pamogk-som-expressed-pw_id={pw_id}.gpickle'
        self.get_cnv_pw_path = lambda \
                pw_id, cnv_type: self.kernel_dir / f'pamogk-cnv-{cnv_type}-pw_id={pw_id}.gpickle'

    @timeit
    def read_som_data(self):
        """
        Returns
        -------
        mapping of patient to mutations by entrez ids
        """
        # Real Data #
        # process RNA-seq expression data
        patients = {}
        with open(config.get_safe_data_file(self.args.som_patient_data)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                ent_id = row['Entrez Gene ID']
                if pat_id not in patients:
                    patients[pat_id] = {ent_id}
                else:
                    patients[pat_id].add(ent_id)

        return collections.OrderedDict(sorted(patients.items()))

    @timeit
    def read_cnv_data(self):
        """
        Returns
        -------
        mapping of patient to mutations by entrez ids
        """
        # Real Data #
        # process Copy Number Variation data
        patients = {}
        ENTREZ_COL_NAME = 'Entrez Gene ID'
        with open(config.get_safe_data_file(self.args.cnv_patient_data)) as csvfile:
            for row in csv.DictReader(csvfile):
                ent_id = row[ENTREZ_COL_NAME]
                for k, v in row.items():
                    if k == ENTREZ_COL_NAME or v == '0':
                        continue
                    pat_id = k
                    if pat_id not in patients:
                        patients[pat_id] = {'loss': [], 'gain': []}
                    if v == '1':
                        patients[pat_id]['loss'].append(ent_id)
                    else:
                        patients[pat_id]['gain'].append(ent_id)

        # map gene ids to uniprot
        # self.preprocess_cnv_patient_data(patients)
        # return as key ordered dict
        return collections.OrderedDict(sorted(patients.items()))

    @timeit
    def find_intersection_patients(self, som_pat, cnv_pat):
        som_pat_list = simplify_pat_ids(som_pat.keys())
        cnv_pat_list = simplify_pat_ids(cnv_pat.keys())

        intersection_list = list(set(som_pat_list).intersection(cnv_pat_list))
        intersection_list.sort()
        save_csv(self.data_dir / 'patients.csv', [[pid] for pid in intersection_list])

        def clean_mapping(s, intersect):
            ns = {}
            for k, v in s.items():
                if k in intersect:
                    ns[k] = v
            return ns

        som_pat = clean_mapping(som_pat, intersection_list)
        cnv_pat = clean_mapping(cnv_pat, intersection_list)

        return som_pat, cnv_pat

    @timeit
    def preprocess_som_patient_data(self, patients):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        res = []
        for pat_id, ent_ids in patients.items():
            # uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
            # if there are any matches map them
            res.append({
                'pat_id': pat_id,
                'mutated_nodes': uni_ids,
            })

        return res

    @timeit
    def preprocess_cnv_patient_data(self, patients):
        # get the dictionary of gene id mappers
        uni2ent, ent2uni = uniprot_mapper.json_to_dict()

        for patient in patients.values():
            for cnv_type, ent_ids in patient.items():
                patient[cnv_type] = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]

    @timeit
    def read_pathways(self):
        # get all pathways
        return cx_pw.read_pathways()

    def som_pathways_save_valid(self, all_pw_map):
        return np.all([self.get_som_pw_path(pw_id).exists() for pw_id in all_pw_map])

    def cnv_pathways_save_valid(self, all_pw_map, cnv_type):
        return np.all([self.get_cnv_pw_path(pw_id, cnv_type).exists() for pw_id in all_pw_map])

    @timeit
    def restore_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_som_pw_path(pw_id)
            logr(f'Loading somatic mutation data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def restore_cnv_pathways(self, all_pw_map, cnv_type):
        num_pw = len(all_pw_map)
        res_pw_map = collections.OrderedDict()
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = self.get_cnv_pw_path(pw_id, cnv_type)
            logr(f'Loading cnv {cnv_type} data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            res_pw_map[pw_id] = nx.read_gpickle(path)
        log()
        return res_pw_map

    @timeit
    def save_som_pathways(self, all_pw_map):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_som_pw_path(pw_id)
            logr(f'Saving somatic mutation data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            nx.write_gpickle(pw, path)
        log()

    @timeit
    def save_cnv_pathways(self, all_pw_map, cnv_type):
        num_pw = len(all_pw_map)
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = self.get_cnv_pw_path(pw_id, cnv_type)
            logr(f'Saving cnv {cnv_type} data {ind + 1:3}/{num_pw} pw_id={pw_id}')
            nx.write_gpickle(pw, path)
        log()

    def label_som_patient_genes(self, all_pw_map, patients):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        patients: :obj:`list`
            list of patients with mutation mappings
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
            logr(f'Checking patient for somatic mutation {ind + 1:4}/{num_pat} pid={pid}')
            label_mapper.mark_label_on_pathways('som', pid, all_pw_map, genes, self.label)
        log()

        self.save_som_pathways(all_pw_map)
        return all_pw_map

    def label_cnv_patient_genes(self, all_pw_map, patients, cnv_type):
        """Labels all patients with matching level of expression

        Parameters
        ----------
        all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
            a dictionary of all pathways we are using
        patients: :obj:`dict`
            mapping of patients with mutation mappings
        """
        # check if we already stored all over/under expression pathway data if so restore them
        if self.cnv_pathways_save_valid(all_pw_map, cnv_type):
            return self.restore_cnv_pathways(all_pw_map, cnv_type)

        num_pat = len(patients)
        # if there are missing ones calculate all of them
        log(f'CNV {cnv_type} patient pathway labeling')
        ind = 1
        for pid, patient in patients.items():
            genes = np.array([patient[cnv_type]])  # get uniprot gene ids from indices
            logr(f'Checking patient for cnv {cnv_type} {ind:4}/{num_pat} pid={pid}')
            label_mapper.mark_label_on_pathways(f'cnv_{cnv_type}', pid, all_pw_map, genes, self.label)
            ind += 1
        log()

        self.save_cnv_pathways(all_pw_map, cnv_type)
        return all_pw_map

    @timeit
    def create_seq_kernels(self, all_pw_map, pat_ids, kms_file_name):
        # experiment variables
        num_pat = pat_ids.shape[0]
        num_pw = len(all_pw_map)
        kms_path = self.kernel_dir / f'{kms_file_name}.npz'
        if kms_path.exists(): return np_load_data(kms_path, key='kms')
        # calculate kernel matrices for over expressed genes
        over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            over_exp_kms[ind] = kernel(pat_ids, pw, label_key=f'label-oe-{self.label}', alpha=self.smoothing_alpha,
                                       normalization=self.kernel_normalization)
            logr(f'Calculating oe pathway kernel={kms_file_name} {ind + 1:4}/{num_pw} pw_id={pw_id}')
        log()

        # calculate kernel matrices for under expressed genes
        under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            under_exp_kms[ind] = kernel(pat_ids, pw, label_key=f'label-ue-{self.label}', alpha=self.smoothing_alpha,
                                        normalization=self.kernel_normalization)
            logr(f'Calculating ue pathway kernel={kms_file_name} {ind + 1:4}/{num_pw} pw_id={pw_id}')
        log()

        kms = np.vstack([over_exp_kms, under_exp_kms])  # stack all kernels
        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def create_cnv_kernels(self, all_pw_map, patients, cnv_type):
        # experiment variables
        num_pat = len(patients)
        num_pw = len(all_pw_map)
        kms_path = self.kernel_dir / f'cnv-kms-{cnv_type}.npz'
        if kms_path.exists(): return np_load_data(kms_path, key='kms')
        # calculate kernel matrices for over expressed genes
        kms = np.zeros((num_pw, num_pat, num_pat))
        pat_ids = np.array(list(patients.keys()))
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = kernel(pat_ids, pw, label_key=f'label-cnv_{cnv_type}', alpha=self.smoothing_alpha,
                              normalization=self.kernel_normalization)
            logr(f'Calculating cnv {cnv_type} pathway kernel {ind + 1:4}/{num_pw} pw_id={pw_id}')
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def create_som_kernels(self, all_pw_map, patients):
        # experiment variables
        num_pat = len(patients)
        num_pw = len(all_pw_map)
        kms_path = self.kernel_dir / 'som-kms.npz'
        if kms_path.exists(): return np_load_data(kms_path, key='kms')
        # calculate kernel matrices for over expressed genes
        kms = np.zeros((num_pw, num_pat, num_pat))
        pat_ids = np.array([pat['pat_id'] for pat in patients])
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):  # for each pathway
            kms[ind] = kernel(pat_ids, pw, label_key='label-som', alpha=self.smoothing_alpha,
                              normalization=self.kernel_normalization)
            logr(f'Calculating som mut pathway kernel {ind + 1:4}/{num_pat} pw_id={pw_id}')
        log()

        np.savez_compressed(kms_path, kms=kms)  # save kernels

        return kms

    @timeit
    def cluster_multiple(self, cluster_sizes, kernel_path):
        if type(cluster_sizes) != list:
            raise TypeError('Cluster sizes must be a list')
        for k in cluster_sizes:
            log(f'Running clustering for k={k}')
            self.cluster(k, kernel_path)

    @timeit
    def cluster(self, n_clusters, kernel_path):
        typ=''
        log('Clustering with smspk-cont')
        # return
        typ_c = ''
        if typ != '':
            typ_c = typ + '_'
        # Cluster using Mkkm-MR
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
        eng.addpath(str(self.data_dir))
        # sending input to the function
        # eng.smspk_clustering_drop_fnc(str(self.exp_data_dir), cluster, drop_percent, typ)
        eng.smspk_clustering_fnc(str(kernel_path), str(self.result_dir), n_clusters, self.drop_percent, typ)
        log('MKKM-MR and K-Means done.')

    @timeit
    def run(self):
        # Somatic mutation data
        som_patients = self.read_som_data()

        # Copy Number Variation Data
        cnv_patients = self.read_cnv_data()

        # Find intersect
        # we don't really use somatic data but return it in case we add it, intersect generated patient.csv that
        # is used by label mapper
        intersect = self.find_intersection_patients(som_patients, cnv_patients)
        som_patients, cnv_patients = intersect

        # CNV loss/gain Data
        all_cnv_pw_map = self.read_pathways()
        labeled_loss_cnv_pw_map = self.label_cnv_patient_genes(all_cnv_pw_map, cnv_patients, 'loss')
        cnv_kernels_loss = self.create_cnv_kernels(labeled_loss_cnv_pw_map, cnv_patients, 'loss')

        all_cnv_pw_map = self.read_pathways()
        labeled_gain_cnv_pw_map = self.label_cnv_patient_genes(all_cnv_pw_map, cnv_patients, 'gain')
        cnv_kernels_gain = self.create_cnv_kernels(labeled_gain_cnv_pw_map, cnv_patients, 'gain')

        kernels = np.concatenate((cnv_kernels_loss, cnv_kernels_gain))

        total = kernels.shape[1] * kernels.shape[2]
        limit = (self.drop_percent * total) / 100.0
        valid_kernels = kernels[np.count_nonzero(kernels, axis=(1, 2)) >= limit]

        log(f'kernel_count={kernels.shape[0]} valid_kernel_count={valid_kernels.shape[0]}')
        kernel_path = self.data_dir / f'kernels-dropped={self.drop_percent}.mat'
        scipy.io.savemat(kernel_path, { 'kernels': valid_kernels})

        cluster_sizes = [2, 3, 4, 5]
        self.cluster_multiple(cluster_sizes, kernel_path)

        self.label_analyzer = LabelAnalysis(results_dir=self.result_dir, methods=['mkkm', 'kmeans'],
                                            cluster_sizes=cluster_sizes, log2_lambdas=self.log2_lambdas)
        self.label_analyzer.run()


def create_experiment(*nargs):
    global args

    if __name__ == '__main__':  # if running directly use command line arguments
        args = parser.parse_args()
    else:  # otherwise use user given arguments
        args = parser.parse_args(nargs)

    print_args(args)

    return Experiment1(args)


if __name__ == '__main__':
    create_experiment().run()
