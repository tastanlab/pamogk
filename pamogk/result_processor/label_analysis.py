#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Label Evaluation
"""
import json

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from pamogk import config
from pamogk.lib.sutils import *
from pamogk.result_processor.latex_generator import pandas_to_latex_table

example_exp_path = config.DATA_DIR / 'pamogk_kirc' / 'pamogk_exp-Experiment1' / \
                   'label=th196-smoothing_alpha=0.9-kr_norm=True'
default_clinical_data_path = config.DATA_DIR / 'kirc_data/kirc_clinical_data.csv'

# create the parser with some defaults
parser = argparse.ArgumentParser(description='Evaluate labels')
parser.add_argument('--clinical-data', '-c', metavar='file-path', dest='clinical_data_path', type=str2path,
                    help='Clinical Data', default=default_clinical_data_path)
parser.add_argument('--show-conf-int', '-ci', dest='show_ci', action='store_true', help='Show Confidence Intervals')
parser.add_argument('--exp-data-dir', '-e', metavar='directory', dest='exp_data_dir', type=str2path,
                    help='Experiment Data Dir', default=example_exp_path)

args = None


class LabelAnalysis(object):
    def __init__(self, exp_data_dir, show_ci=False, clinical_data_path=default_clinical_data_path, methods=None,
                 cluster_sizes=None, log2_lambdas=None):
        if log2_lambdas is None:
            log2_lambdas = list(range(-15, 16, 3))
        if cluster_sizes is None:
            cluster_sizes = list(range(2, 6))
        if methods is None:
            methods = ['mkkm', 'kmeans']
        print_args(locals())

        self.exp_data_dir = Path(exp_data_dir)
        self.show_ci = show_ci
        self.res_dir = self.exp_data_dir / 'results'
        self.fig_dir = self.exp_data_dir / 'figures'
        self.exports_path = self.fig_dir / 'exports.json'
        self.log2_lambdas = log2_lambdas
        self.cluster_sizes = cluster_sizes
        self.methods = methods

        safe_create_dir(self.fig_dir)
        # load clinical and experiment data
        clinical_data = self.read_clinical_data(clinical_data_path)
        exp_pat_ids = self.load_exp_patients()
        # this array holds true for patients that have data for the patient at this index
        self.exp_pats_with_data = np.array([pat_id in clinical_data for pat_id in exp_pat_ids])
        # these are clinical data for patients with data
        pat_ids = exp_pat_ids[self.exp_pats_with_data]
        self.clinical_days = np.array([clinical_data[pat_id]['days'] for pat_id in pat_ids]).astype(int)
        self.clinical_status = np.array([clinical_data[pat_id]['status'] for pat_id in pat_ids]).astype(int)
        self.exported_files = {}

    @timeit
    def load_exp_patients(self):
        """
        Loads list of intersecting patients from exp data dir
        """
        with open(self.exp_data_dir / 'patients.csv') as f:
            reader = csv.reader(f)
            return np.array([row[0] for row in reader])

    @timeit
    def read_label_data(self, filename):
        """
        Reads the label data and patient_id data and create array of dict with elements pat_id and label
        Parameters
        ----------
        filename

        Returns
        -------

        """
        return np_load_data(self.res_dir / f'{filename}.npz', key='labels')

    @timeit
    def read_clinical_data(self, file_loc):
        """
        Load read clinical data
        Parameters
        ----------
        file_loc

        Returns
        -------

        """
        clinical_data = {}
        with open(file_loc) as f:
            reader = csv.DictReader(f)
            for row in reader:
                clinical_data[row['Patient ID']] = dict(status=row['Status'], days=row['Days'])
        return clinical_data

    def calc_logrank_p_value(self, labels):
        """
        Calculates p value for given labels using clinical data
        Parameters
        ----------
        labels: iterable labels, must have same size as self.exp_pats_with_data

        Returns
        -------

        """
        return multivariate_logrank_test(self.clinical_days, labels, self.clinical_status).p_value

    @timeit
    def km_analysis(self, labels, out_filename):
        """
        Parameters
        ----------
        labels:
            Array patient dicts, which includes pat_id, label, status, days
        out_filename:
            Name of file to save Kaplan Meier Graph

        Return
        -----------
        [0]:
            Multivariate log-rank test p values
        [1]:
            One-vs-All result of log-rank test (array of p values)
        """
        plt.clf()

        ax = plt.subplot(111)
        ax.set_xlabel('Time (day)')
        ax.set_ylabel('Survival Probability')
        for label in np.unique(labels):  # get unique cluster labels
            cluster_pat_ind = labels == label  # store indices of patients in this cluster
            days = self.clinical_days[cluster_pat_ind]
            status = self.clinical_status[cluster_pat_ind]
            kmf = KaplanMeierFitter()
            kmf.fit(days, status, label=f'size={np.sum(cluster_pat_ind)}')
            # save to be used by others
            out_path = self.fig_dir / f'{out_filename}-label={label}.npz'
            np_save_npz(out_path, kmf_prob_func=kmf.survival_function_.values,
                        kmf_conf_int=kmf.confidence_interval_.values, kmf_timeline=kmf.timeline)
            self.add_exported_filepath(out_path)
            kmf.plot(ax=ax, ci_show=self.show_ci)

        out_path = self.fig_dir / f'{out_filename}.png'
        plt.savefig(out_path)
        self.add_exported_filepath(out_path)

        # all
        all_p = self.calc_logrank_p_value(labels)
        # one vs all
        vs_p = [f"{self.calc_logrank_p_value(labels == lb):.2e}" for lb in np.unique(labels)]

        return all_p, vs_p

    # Process one file
    def process_label_file(self, filename):
        """
        Takes label from label_file and saves figure to out_file location
        Returns the multivariate p result and one-vs-all p result.

        Parameters
        ----------
        filename: str
            Filename that shows the path for the kernel

        Return
        -----------
        [0]:
            Multivariate log-rank test p values
        [1]:
            One-vs-All result of log-rank test (array of p values)
        """
        labels = np_load_data(self.res_dir / f'{filename}.npz', key='labels')
        # filter out patients that have no data
        labels = labels[self.exp_pats_with_data]

        p, vs_p = self.km_analysis(labels, filename)

        return p, vs_p

    def add_exported_filepath(self, path):
        path = Path(path).absolute()
        if not path.exists():
            ValueError('Tried to add non existing path to exports')
        if path.suffix not in self.exported_files:
            self.exported_files[path.suffix] = []
        self.exported_files[path.suffix].append(str(path))

    def df_to_csv(self, df, name):
        out_path = (self.fig_dir / name).with_suffix('.csv')
        df.to_csv(out_path)
        self.add_exported_filepath(out_path)

    def df_to_latex_table(self, df, name, row_name='l', col_name='k'):
        out_path = (self.fig_dir / name).with_suffix('.tex')
        pandas_to_latex_table(df, row_name, col_name, out_path)
        self.add_exported_filepath(out_path)

    @timeit
    def run(self):
        log2_lmbds_str = [str(lmbd) for lmbd in self.log2_lambdas]
        result_inds = log2_lmbds_str.copy()

        for method in self.methods:
            result_inds.append(method)

        result_df = pd.DataFrame(columns=self.cluster_sizes, index=result_inds)

        for method in self.methods:
            for label in self.cluster_sizes:
                # compare dataframe
                label_i = int(label)
                comp_df = pd.DataFrame(columns=list(range(1, label_i + 1)), index=log2_lmbds_str)
                if method == 'mkkm':
                    for (log2_lmbd, log2_lmbd_str) in zip(self.log2_lambdas, log2_lmbds_str):
                        res, vs_res = self.process_label_file(f'pamogk-{method}-k={label}-log2_lambda={log2_lmbd}')
                        result_df.loc[log2_lmbd_str][label] = f'{res:.2e}'
                        comp_df.loc[log2_lmbd_str] = vs_res

                    self.df_to_latex_table(comp_df, f'latex_table_1v_{label_i - 1}.tex')
                    self.df_to_csv(comp_df, f'results_1v_{label_i - 1}')
                else:
                    res, vs_res = self.process_label_file(f'pamogk-{method}-k={label}')
                    result_df.loc[method][label] = f'{res:.2e}'

        self.df_to_latex_table(result_df, 'latex_table')
        self.df_to_csv(result_df, 'results')
        with open(self.exports_path, 'w') as f:
            json.dump(self.exported_files, f)
        log(f'Finished label analysis with exported files on path={self.exports_path}')


if __name__ == '__main__':
    # if running directly use command line arguments
    args = parser.parse_args()
    LabelAnalysis(**vars(args)).run()
