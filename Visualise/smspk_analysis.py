#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import csv
from lib.sutils import *
import argparse
from data_processor import rnaseq_processor as rp
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--rnaSeq-data', '-r', metavar='file-path', dest='rnaseq_data', type=str, help='RnaSeq Data', default='../data/kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')
parser.add_argument('--clinical-data', '-c', metavar='file-path', dest='clinical_data', type=str, help='Clinical Data', default='../data/kirc_data/kirc_clinical_data.csv')
parser.add_argument('--label-data', '-l', metavar='file-path', dest='label_data', type=str, help='Label Data', default='../experiments/label_data/smspk-2lab-5it-lmkkmeans.txt')

args = parser.parse_args()
log('Running args:', args)


class Experiment1(object):
    def __init__(self):
        print()

    @timeit
    def read_data(self):
        ### Real Data ###
        # process RNA-seq expression data
        gene_exp, gene_name_map = rp.process(args.rnaseq_data)

        # convert entrez gene id to uniprot id
        pat_ids = gene_exp.columns.values # patient TCGA ids

        for i in range(len(pat_ids)):
            pat_ids[i] = "-".join(pat_ids[i].split("-")[0:3])
        return pat_ids

    @timeit
    def labelize(self,pat_ids):
        patients=[]
        labels = np.loadtxt(args.label_data)
        for i in range(len(pat_ids)):
            patients.append({"pat_id":pat_ids[i],"label":labels[i]})

        return patients

    @timeit
    def read_clinical_data(self,fileLoc):
        ### Real Data ###
        # process RNA-seq expression data
        patients = {}
        with open(fileLoc) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pat_id = row['Patient ID']
                status = row['Status']
                days = row['Days']
                patients[pat_id] = [status,days]

        return patients

    @timeit
    def findIntersection(self,patients, clinicalData):
        ### Real Data ###
        # Eliminate patients not in ClinicalData, insert status and days
        delIndex = []
        for i in range(len(patients)):
            if patients[i]["pat_id"] in clinicalData:
                patients[i]["status"] = clinicalData[patients[i]["pat_id"]][0]
                patients[i]["days"] = clinicalData[patients[i]["pat_id"]][1]
            else:
                print("Patient "+str(patients[i]["pat_id"])+" will be deleted")
                delIndex.append(i)

        return np.delete(patients,delIndex)

    @timeit
    def KMAnalysis(self,patients):

        labels = np.array([p["label"] for p in patients])
        days = np.array([p["days"] for p in patients])
        status = np.array([p["status"] for p in patients])

        clusterSize = len(np.unique(labels))
        ax = plt.subplot(111)
        for i in range(clusterSize):
            indexes = [index for index in range(len(patients)) if labels[index]==i]
            daysTmp = days[indexes]
            statusTmp = status[indexes]
            kmf = KaplanMeierFitter()
            kmf.fit(daysTmp, statusTmp, label="Cluster"+str(i))
            kmf.plot(ax=ax)
        plt.show()
        df = pd.DataFrame({
            'durations': days.astype(float),
            'groups': labels,
            'events': status.astype(int),
        })

        results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])

        return results.p_value


exp = Experiment1()

patient_map = exp.read_data()

patients = exp.labelize(patient_map)

clinicalData = exp.read_clinical_data(args.clinical_data)

patients = exp.findIntersection(patients, clinicalData)

p = exp.KMAnalysis(patients)

print("P value: ",p)
