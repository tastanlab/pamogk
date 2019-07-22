#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper
import csv
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=str, help='pathway ID list', default='../data/kirc_data/kirc_somatic_mutation_data.csv')
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable Debug Mode')
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float, help='Node2Vec feature space size', default=128)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed', default=False)
parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments', default=1000)
parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float, help='Surviving patient percentage in range [0, 1]', default=0.9)
parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]', default=0.4)

args = parser.parse_args()


def read_data():
    ### Real Data ###
    # process RNA-seq expression data
    patients = {}
    with open(args.patient_data) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pat_id = row['Patient ID']
            ent_id = row['Entrez Gene ID']
            if pat_id not in patients: patients[pat_id] = set([ent_id])
            else: patients[pat_id].add(ent_id)

    return patients

def preprocess_patient_data(patients):
    # get the dictionary of gene id mappers
    uni2ent, ent2uni = uniprot_mapper.json_to_dict()

    res = []
    num_empty = 0
    for pat_id, ent_ids in patients.items():
        # uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
        uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
        # if there are any matches map them
        if len(uni_ids) > 0: res.append({
            'pat_id': pat_id,
            'mutated_nodes': uni_ids,
        })
        else: num_empty += 1


    return res

def read_pathways():
    # get all pathways
    return cx_pw.read_pathways()

def get_nodes(all_pw_map):
    res = {}
    num_pw = len(all_pw_map)
    for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
        gene_list = []
        for node in pw._node:
            gene_list += pw._node[node]["uniprot-ids"]
        res[pw_id] = gene_list

    return res

def get_stats_from_pathways(patients, pathway_gene_map):
    no_of_pathways = len(pathway_gene_map)
    patient_not_count_in_pathways = []
    for patient in patients:
        hit = 0
        not_hit = 0
        for pathway in pathway_gene_map.keys():
            if len(np.intersect1d(pathway_gene_map[pathway], patient["mutated_nodes"])) == 0:
                not_hit += 1
            else:
                hit += 1
        patient["hit"]=hit
        patient["not_hit"]=not_hit
        patient_not_count_in_pathways.append(not_hit)

    return patient_not_count_in_pathways


patient_map = read_data()

patients = preprocess_patient_data(patient_map)

all_pw_map = read_pathways()

pathway_gene_map = get_nodes(all_pw_map)

patient_pathway_stats = get_stats_from_pathways(patients, pathway_gene_map)

plt.hist(patient_pathway_stats)  # arguments are passed to np.histogram
plt.title("Histogram with 'auto' bins")
plt.show()

print()

