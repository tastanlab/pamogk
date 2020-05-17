#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Label Evaluation
Arguments take 3 files:
    patient_data: list of patients (csv file 1 row)
    clinical_data: csv clinical data file.
    label_file: file containing labels which can be read with np.loadtxt() function at once.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from pamogk import config
from . import report_creator as rc
from pamogk.lib.sutils import *

KERNELS = ['pamogk-all']
METHODS = ['mkkm', 'kmeans']
LABELS = ['2', '3', '4', '5']
LAMBDAS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']

parser = argparse.ArgumentParser(description='Evaluate labels')
parser.add_argument('--patient-data', '-p', metavar='file-path', dest='patient_data', type=Path, help='Patient Ids',
                    default=config.DATA_DIR / 'kirc_data/kirc_intersect.csv')
parser.add_argument('--clinical-data', '-c', metavar='file-path', dest='clinical_data', type=Path, help='Clinical Data',
                    default=config.DATA_DIR / 'kirc_data/kirc_clinical_data.csv')
parser.add_argument('--label-file', '-l', metavar='file-path', dest='label_file', type=Path, help='Label Data',
                    default=config.DATA_DIR / 'pamogk_all/Experiment1-label=1-smoothing_alpha=0.9-norm=True/labels/')

args = parser.parse_args()
print_args(args)


class LabelAnalysis(object):
    def __init__(self):
        print()

    # Reads the label data and patient_id data and create array of dict with elements pat_id and label
    @timeit
    def read_data(self, label_loc):
        patients = []
        labels = np.loadtxt(label_loc)
        with open(args.patient_data, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            header = next(reader)
            for idx, pt_id in enumerate(header):
                patients.append({'pat_id': pt_id, 'label': labels[idx]})

        return patients

    @timeit
    def read_clinical_data(self, file_loc):
        # Real Data #
        # process clinical data
        patients = {}
        with open(file_loc) as f:
            reader = csv.DictReader(f)
            for row in reader:
                pat_id = row['Patient ID']
                status = row['Status']
                days = row['Days']
                patients[pat_id] = [status, days]

        return patients

    @timeit
    def findIntersection(self, patients, clinical_data):
        # Real Data #
        # Eliminate patients not in ClinicalData, insert status and days
        del_index = []
        for i in range(len(patients)):
            if patients[i]['pat_id'] in clinical_data.keys():
                patients[i]['status'] = clinical_data[patients[i]['pat_id']][0]
                patients[i]['days'] = clinical_data[patients[i]['pat_id']][1]
            else:
                print('Patient', patients[i]['pat_id'], 'will be deleted')
                del_index.append(i)

        return np.delete(patients, del_index)

    @timeit
    def one_vs_all(self, data_df):
        labels = list(data_df['groups'])
        min_label = min(labels)
        max_label = max(labels)
        vs_p = []
        for cur_label in range(min_label, max_label + 1):
            new_labels = []
            for i, label in enumerate(labels):
                if label == cur_label:
                    new_labels.append(cur_label)
                else:
                    new_labels.append(cur_label + 1)
            results = multivariate_logrank_test(data_df['durations'], new_labels, data_df['events'])
            p = f'{results.p_value:.2e}'
            vs_p.append(p)
        return vs_p

    @timeit
    def km_analysis(self, patients, out_file):
        """
        Parameters
        ----------
        patients:
            Array patient dicts, which includes pat_id, label, status, days
        out_file:
            Location to save Kaplan Meier Graph

        Return
        -----------
        [0]:
            Multivariate log-rank test p values
        [1]:
            One-vs-All result of log-rank test (array of p values)
        """
        plt.clf()
        labels = np.array([p['label'] for p in patients]).astype(int)
        days = np.array([p['days'] for p in patients]).astype(int)
        status = np.array([p['status'] for p in patients]).astype(int)

        ax = plt.subplot(111)
        for i in range(min(labels), max(labels) + 1):
            indexes = [index for index in range(len(patients)) if labels[index] == i]
            days_tmp = days[indexes]
            status_tmp = status[indexes]
            kmf = KaplanMeierFitter()
            kmf.fit(days_tmp, status_tmp, label=str(len(indexes)))
            kmf.plot(ax=ax, ci_show=False)
            ax.set_xlabel('Süre (gün)')
            ax.set_ylabel('Hayatta Kalma Olasılığı')

        plt.savefig(out_file)
        df = pd.DataFrame({
            'durations': days,
            'groups': labels,
            'events': status,
        })

        results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
        vs_p = self.one_vs_all(df)

        return results.p_value, vs_p


# Process one file
def process_one_file(label_file='', out_file='--'):
    """
    Takes label from label_file and saves figure to out_file location
    Returns the multivariate p result and one-vs-all p result.

    Parameters
    ----------
    label_file:
        Label file location, if none given it gets from argument parser
    out_file:
        Output graph location, if none given it saves to label_file folder as no-name.png

    Return
    -----------
    [0]:
        Multivariate log-rank test p values
    [1]:
        One-vs-All result of log-rank test (array of p values)
    """

    if label_file == '':
        label_file = args.label_file
    label_file = Path(label_file)
    if out_file == '--':
        out_file = label_file.parent / 'no-name.png'

    exp = LabelAnalysis()

    patients = exp.read_data(label_file)

    clinical_data = exp.read_clinical_data(args.clinical_data)

    patients = exp.findIntersection(patients, clinical_data)

    p, vs_p = exp.km_analysis(patients, out_file)

    return p, vs_p


def main():
    lambda_values = [2 ** (-18 + 3 * int(lmb)) for lmb in LAMBDAS]
    lambda_values_formatted = [f'{lmb:.2e}' for lmb in lambda_values]
    lambda_values_formatted_for_vs = lambda_values_formatted.copy()
    if 'kmeans' in METHODS:
        lambda_values_formatted.append('kmeans')
    if 'lmkkmeans' in METHODS:
        lambda_values_formatted.append('lmkkmeans')
    result_df = pd.DataFrame(columns=LABELS, index=lambda_values_formatted)

    fig_dir = args.label_file / 'figures'
    safe_create_dir(fig_dir)

    for kernel in KERNELS:
        for method in METHODS:
            for label in LABELS:
                vs_df = pd.DataFrame(columns=list(range(1, int(label) + 1)), index=lambda_values_formatted_for_vs)
                if method == 'mkkm':
                    method_name = 'rmmk'
                    for lamd in LAMBDAS:
                        in_file = args.label_file / f'{kernel}-{method_name}-{label}lab-lambda{lamd}'
                        out_file = fig_dir / f'{kernel}-{method}-{label}lab-lambda{lamd}.png'
                        res, vs_res = process_one_file(in_file, out_file)
                        lambda_value = 2 ** (-18 + 3 * int(lamd))
                        lambda_value_formatted = f'{lambda_value:.2e}'
                        result_df.loc[lambda_value_formatted][label] = f'{res:.2e}'
                        vs_df.loc[lambda_value_formatted] = vs_res

                    out_file = args.label_file / f'latex_table_1v_{int(label) - 1}.txt'
                    rc.pandas_to_latex_table(vs_df, 'l', 'k', out_file)
                    out_file = args.label_file / f'results_1v_{int(label) - 1}.csv'
                    vs_df.to_csv(out_file)
                else:
                    in_file = args.label_file / f'{kernel}-{method}-{label}lab'
                    out_file = fig_dir / f'{kernel}-{method}-{label}lab.png'
                    res, vs_res = process_one_file(in_file, out_file)
                    result_df.loc[method][label] = f'{res:.2e}'

    rc.pandas_to_latex_table(result_df, 'l', 'k', args.label_file / 'latex_table.txt')
    result_df.to_csv(args.label_file / 'results.csv')


if __name__ == '__main__':
    main()
