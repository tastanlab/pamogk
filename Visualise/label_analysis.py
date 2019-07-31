#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
sys.path.append('..')
import matplotlib.pyplot as plt
import pandas as pd
import csv
from lib.sutils import *
import argparse
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from collections import OrderedDict
import Visualise.report_creator as rc

'''
Label Evaluation
Arguments take 3 files:
    patient_data: list of patients (csv file 1 row)
    clinical_data: csv clinical data file.
    label_file: file containing labels which can be read with np.loadtxt() function at once.
'''

parser = argparse.ArgumentParser(description='Evaluate labels')
parser.add_argument('--patient-data', '-p', metavar='file-path', dest='patient_data', type=str, help='Patient Ids',
                    default='../data/kirc_data/kirc_intersect.csv')
parser.add_argument('--clinical-data', '-c', metavar='file-path', dest='clinical_data', type=str, help='Clinical Data',
                    default='../data/kirc_data/kirc_clinical_data.csv')
parser.add_argument('--label-file', '-l', metavar='file-path', dest='label_file', type=str, help='Label Data',
                    default='../data/smspk_all/Experiment1-label=1-smoothing_alpha=0.9-norm=True/labels/')

args = parser.parse_args()
log('Running args:', args)


class Label_Analysis(object):
    def __init__(self):
        print()

    # Reads the label data and patient_id data and create array of dict with elements pat_id and label
    @timeit
    def read_data(self, label_loc):
        patients = []
        labels = np.loadtxt(label_loc)
        with open(args.patient_data, "r") as f:
            reader = csv.reader(f, delimiter=",")
            row1 = next(reader)
            for idx, pt_id in enumerate(row1):
                patients.append({"pat_id": pt_id, "label": labels[idx]})

        return patients

    @timeit
    def read_clinical_data(self, file_loc):
        ### Real Data ###
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
        ### Real Data ###
        # Eliminate patients not in ClinicalData, insert status and days
        delIndex = []
        for i in range(len(patients)):
            if patients[i]["pat_id"] in clinical_data.keys():
                patients[i]["status"] = clinical_data[patients[i]["pat_id"]][0]
                patients[i]["days"] = clinical_data[patients[i]["pat_id"]][1]
            else:
                print("Patient " + str(patients[i]["pat_id"]) + " will be deleted")
                delIndex.append(i)

        return np.delete(patients, delIndex)

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
            p = '{:.2e}'.format(results.p_value)
            vs_p.append(p)
        return vs_p

    @timeit
    def km_analysis(self, patients, out_file):
        '''
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
        '''
        plt.clf()
        labels = np.array([p["label"] for p in patients]).astype(int)
        days = np.array([p["days"] for p in patients]).astype(int)
        status = np.array([p["status"] for p in patients]).astype(int)

        ax = plt.subplot(111)
        for i in range(min(labels), max(labels) + 1):
            indexes = [index for index in range(len(patients)) if labels[index] == i]
            daysTmp = days[indexes]
            statusTmp = status[indexes]
            kmf = KaplanMeierFitter()
            kmf.fit(daysTmp, statusTmp, label=str(len(indexes)))
            kmf.plot(ax=ax, ci_show=False)
            ax.set_xlabel("Süre (gün)")
            ax.set_ylabel("Hayatta Kalma Olasılığı")

        plt.savefig(out_file)
        df = pd.DataFrame({
            'durations': days,
            'groups': labels,
            'events': status,
        })

        results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
        vs_p = self.one_vs_all(df)

        return results.p_value, vs_p

#Process one file
def process_one_file(label_file='', out_file='--'):
    '''
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
    '''

    if label_file == '':
        label_file = args.label_file
    if out_file == '--':
        out_file = os.path.join(os.path.dirname(label_file), 'no-name.png')

    exp = Label_Analysis()

    patients = exp.read_data(label_file)

    clinical_data = exp.read_clinical_data(args.clinical_data)

    patients = exp.findIntersection(patients, clinical_data)

    p, vs_p = exp.km_analysis(patients, out_file)

    return p, vs_p


def main():
    KERNELS = ["smspk-all"]
    METHODS = ["mkkm", "kmeans"]
    LABELS = ["2", "3", "4", "5"]
    LAMBDAS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]

    lambda_values = [2 ** (-18 + 3 * int(l)) for l in LAMBDAS]
    lambda_values_formatted = ['{:.2e}'.format(l) for l in lambda_values]
    lambda_values_formatted_for_vs = lambda_values_formatted.copy()
    if "kmeans" in METHODS:
        lambda_values_formatted.append("kmeans")
    if "lmkkmeans" in METHODS:
        lambda_values_formatted.append("lmkkmeans")
    result_df = pd.DataFrame(columns=LABELS, index=lambda_values_formatted)

    if not os.path.exists(args.label_file + "figures"):
        os.makedirs(args.label_file + "figures")

    for kernel in KERNELS:
        for method in METHODS:
            for label in LABELS:
                vs_df = pd.DataFrame(columns=list(range(1, int(label) + 1)), index=lambda_values_formatted_for_vs)
                if method == "mkkm":
                    method_name = "rmmk"
                    for lambdaa in LAMBDAS:
                        in_file = args.label_file + kernel + "-" + method_name + "-" + label + "lab-lambda" + lambdaa
                        out_file = args.label_file + "figures/" + kernel + "-" + method + "-" + label + "lab-lambda" + lambdaa + ".png"
                        res, vs_res = process_one_file(in_file, out_file)
                        lambda_value = 2 ** (-18 + 3 * int(lambdaa))
                        lambda_value_formatted = '{:.2e}'.format(lambda_value)
                        result_df.loc[lambda_value_formatted][label] = '{:.2e}'.format(res)
                        vs_df.loc[lambda_value_formatted] = vs_res

                    rc.pandas_to_latex_table(vs_df, "l", "k",
                                             args.label_file + "/latex_table_1v_" + str(int(label) - 1) + ".txt")
                    rc.pandas_to_csv_table(vs_df, "l", "k",
                                           args.label_file + "/results_1v_" + str(int(label) - 1) + ".csv")
                else:
                    method_name = method
                    in_file = args.label_file + kernel + "-" + method + "-" + label + "lab"
                    out_file = args.label_file + "figures/" + kernel + '-' + method + "-" + label + "lab.png"
                    res, vs_res = process_one_file(in_file, out_file)
                    result_df.loc[method_name][label] = '{:.2e}'.format(res)

    rc.pandas_to_latex_table(result_df, "l", "k", args.label_file + "/latex_table.txt")
    rc.pandas_to_csv_table(result_df, "l", "k", args.label_file + "/results.csv")


if __name__ == '__main__':
    main()
