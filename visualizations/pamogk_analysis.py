#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

import config
from data_processor import rnaseq_processor as rp
from lib.sutils import *

parser = argparse.ArgumentParser(description='Run PAMOGK algorithms on pathways')
parser.add_argument('--rnaSeq-data', '-r', metavar='file-path', dest='rnaseq_data', type=Path, help='RnaSeq Data',
                    default=config.DATA_DIR / 'kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--clinical-data', '-c', metavar='file-path', dest='clinical_data', type=Path, help='Clinical Data',
                    default=config.DATA_DIR / 'kirc_data/kirc_clinical_data.csv')
parser.add_argument('--label-data', '-l', metavar='file-path', dest='label_data', type=Path, help='Label Data',
                    default=config.ROOT_DIR / 'experiments/label_data/pamogk-2lab-5it-lmkkmeans.txt')

args = parser.parse_args()
print_args(args)


class Experiment1(object):
    def __init__(self):
        print()

    @timeit
    def read_data(self):
        # Real Data #
        # process RNA-seq expression data
        gene_exp, gene_name_map = rp.process(args.rnaseq_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values  # patient TCGA ids

        for i in range(len(pat_ids)):
            pat_ids[i] = '-'.join(pat_ids[i].split('-')[0:3])
        return pat_ids

    @timeit
    def labelize(self, pat_ids):
        patients = []
        labels = np.loadtxt(args.label_data)
        for i in range(len(pat_ids)):
            patients.append({'pat_id': pat_ids[i], 'label': labels[i]})

        return patients

    @timeit
    def read_clinical_data(self, filepath):
        # Real Data #
        # process RNA-seq expression data
        patients = {}
        with open(filepath) as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                pat_id = row['Patient ID']
                status = row['Status']
                days = row['Days']
                patients[pat_id] = [status, days]

        return patients

    @timeit
    def find_intersection(self, patients, clinical_data):
        # Real Data #
        # Eliminate patients not in ClinicalData, insert status and days
        del_index = []
        for i in range(len(patients)):
            if patients[i]['pat_id'] in clinical_data:
                patients[i]['status'] = clinical_data[patients[i]['pat_id']][0]
                patients[i]['days'] = clinical_data[patients[i]['pat_id']][1]
            else:
                print('Patient ' + str(patients[i]['pat_id']) + ' will be deleted')
                del_index.append(i)

        return np.delete(patients, del_index)

    @timeit
    def km_analysis(self, patients):

        labels = np.array([p['label'] for p in patients])
        days = np.array([p['days'] for p in patients])
        status = np.array([p['status'] for p in patients])

        pat_no = []
        colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']
        unique_labels = np.unique(labels)
        # cluster_size = len(np.unique(labels))
        ax = plt.subplot(111, figsize=(5.2, 4.8))
        for idx, i in enumerate(unique_labels):
            indexes = [index for index in range(len(patients)) if labels[index] == i]
            pat_no.append(str(len(indexes)))
            kmf = KaplanMeierFitter()
            kmf.fit(days[indexes], status[indexes], label='Cluster' + str(i))
            kmf.plot(ax=ax, color=colors[idx], ci_show=False, show_censors=True,
                     censor_styles={'ms': 4, 'mew': 0.2, 'marker': '|'})

        ax.set_xlabel("Time (days)", fontsize=14, fontname="Arial")
        ax.set_ylabel("Survival Probability", fontsize=14, fontname="Arial")
        hand = [plt.plot([], ls="-", color=colors[i])[0] for i in range(len(pat_no))]
        labs = pat_no
        plt.legend(handles=hand, labels=labs, frameon=False, title='Number of Patients',
                   loc='upper right')  # prop = {"name":'Arial'}
        plt.yticks(np.arange(0, 1.2, step=0.2))
        plt.ylim(0, 1.1)
        df = pd.DataFrame({
            'durations': days.astype(float),
            'groups': labels,
            'events': status.astype(int),
        })

        results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
        text_p = results.p_value
        plt.text(0.05, 0.08, f'p-value: {text_p:.2e}', transform=ax.transAxes, fontdict={"name": 'Arial'}, fontsize=12,
                 color="black")
        ax.set_xlabel("Time (days)", fontsize=14, fontname="Arial")
        ax.set_ylabel("Survival Probability", fontsize=14, fontname="Arial")
        plt.rcParams['legend.title_fontsize'] = 12
        # plt.savefig(out_file, format="pdf", dpi=600)
        plt.show()

        return results.p_value


def main():
    exp = Experiment1()

    patient_map = exp.read_data()

    patients = exp.labelize(patient_map)

    clinical_data = exp.read_clinical_data(args.clinical_data)

    patients = exp.find_intersection(patients, clinical_data)

    p = exp.km_analysis(patients)

    print('P value: ', p)


if __name__ == '__main__':
    main()
