#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import networkx as nx
import copy
import pdb
from operator import itemgetter
import operator
import matplotlib.pyplot as plt
import smspk
import pandas as pd
import label_mapper
from structural_processor import rnaseq_processor as rp
from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper as um
import time
import json
import os
import config
from lib.sutils import *

### Real Data ###
# process RNA-seq expression data
gene_exp, gene_name_map = rp.process('data/kirc_data/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')

# get the dictionary of gene id mappers
[uni2ent, ent2uni] = um.json_to_dict()

# convert entrez gene id to uniprot id
pat_ids = gene_exp.columns.values # patient TCGA ids
all_ent_ids = gene_exp.index.values # gene entrez ids
GE = gene_exp.values

found_ent_ids = [eid in ent2uni for eid in all_ent_ids]
ent_ids = np.array([eid for eid in all_ent_ids if eid in ent2uni])
uni_ids = np.array([ent2uni[eid] for eid in ent_ids])

print('uni_ids:', len(uni_ids))
print('miss_ent_ids:', len(all_ent_ids) - sum(found_ent_ids))

# prune genes whose uniprot id is not found
GE = GE[found_ent_ids]
# pdb.set_trace()

# get all pathways
all_pw_map = cx_pw.read_pathways()

### Synthetic Data ###
# # synthetic gene expression data
# gene_exp = pd.DataFrame(np.array([[1, -1, -1], [1, 0, 1], [0, 1, -1], [-1, 1, 1], [0, 0, 1]], dtype=int), columns=['pt1', 'pt2', 'pt3'])
# gene_exp['uniprot_id'] = ['A', 'B', 'D', 'F', 'X']
# gene_exp['at_num'] = ['at1', 'at2', 'at3', 'at4', 'at5']
# # gene_exp = pd.DataFrame(np.array([[1, 0, 0, 'A', 'at1'], [1, 0, 1, 'B', 'at2'], [0, 1, 0, 'D', 'at3'], [0, 1, 1, 'F', 'at4'], [0, 0, 1, 'X', 'at5']]), columns=['pt1', 'pt2', 'pt3', 'uniprot_id', 'at_num'])
# gene_exp = gene_exp.set_index('at_num')
#
# # synthetic pathways
# # Create a networkx graph object
# mg = nx.Graph()
# # mg.add_edges_from([(1,2), (1,3), (3,4), (1,5), (3,5), (4,2), (2,3), (3,0)])
# mg.add_edges_from([('A','B'), ('A','C'), ('C','D'), ('A','E'), ('C','E'), ('D','B'), ('B','C'), ('C','X')])
#
# # nx.set_node_attributes(mg, {0:'Protein', 1:'Protein', 2:'Protein', 3:'Calcium', 4:'Protein', 5:'Calcium'}, 'type')
# # nx.set_node_attributes(mg, {0:['uniprot:X'], 1:['uniprot:A'], 2:['uniprot:B'], 3:['uniprot:C'], 4:['uniprot:D'], 5:['uniprot:E']}, 'alias')
# nx.set_node_attributes(mg, {'X':'Protein', 'A':'Protein', 'B':'Protein', 'C':'Calcium', 'D':'Protein', 'E':'Calcium'}, 'type')
# nx.set_node_attributes(mg, {'X':['uniprot:X'], 'A':['uniprot:A'], 'B':['uniprot:B'], 'C':['uniprot:C'], 'D':['uniprot:D'], 'E':['uniprot:E']}, 'alias')
#
# mg2 = nx.Graph()
# # mg2.add_edges_from([(1,3), (1,4), (1,2), (3,2), (4,0)])
# mg2.add_edges_from([('B','D'), ('B','E'), ('B','G'), ('D','G'), ('E','F')])
#
# # nx.set_node_attributes(mg2, {0:'Protein', 1:'Protein', 2:'Calcium', 3:'Protein', 4:'Calcium'}, 'type')
# # nx.set_node_attributes(mg2, {0:['uniprot:F'], 1:['uniprot:B'], 2:['uniprot:G'], 3:['uniprot:D'], 4:['uniprot:E']}, 'alias')
# nx.set_node_attributes(mg2, {'F':'Protein', 'B':'Protein', 'G':'Calcium', 'D':'Protein', 'E':'Calcium'}, 'type')
# nx.set_node_attributes(mg2, {'F':['uniprot:F'], 'B':['uniprot:B'], 'G':['uniprot:G'], 'D':['uniprot:D'], 'E':['uniprot:E']}, 'alias')
# pw_map = {'pw1': mg, 'pw2': mg2}


@timeit
def generate_pat_map(patient_ids, pw_map):
    '''Generates a patient id to pathway to gene mapping

    Parameters
    ----------
    patient_ids: list of int
        List of patient ids showing which columns of gene expression date are
        the patients
    pw_map: map of cx pathway id to networx graph
        Holds the mapping for all pathways belonging to cx pathway graphs
    '''
    return dict((pid, dict((k, {}) for k in pw_map.keys())) for pid in patient_ids)


# pdb.set_trace()
@timeit
def label_patient_genes(all_pw_map, pat_ids, GE, label=1):
    '''Labels all patients with matching level of expression

    Parameters
    ----------
    all_pw_map: :obj:`list` of :obj:`networkx.classes.graph.Graph`
        a dictionary of all pathways we are using
    pat_ids: :obj:`list` of :obj:`str`
        list of patient ids
    GE: :obj:`numpy.ndarray`
        Gene expression data array in shape of genes by patients
    label: int, optional
        label that will be used for marking patients
    '''
    graph_dir = os.path.join(config.data_dir, 'smspk')
    safe_create_dir(graph_dir)
    graph_file = 'smspk-over-under-expressed-label={}'.format(label)
    graph_path = os.path.join(graph_dir, graph_file);

    get_pw_path = lambda pw_id: '{}-pw_id={}.gpickle'.format(graph_path, pw_id)

    num_pw = len(all_pw_map)

    @timeit
    def restore_pathways():
        for ind, pw_id in enumerate(all_pw_map.keys()):
            path = get_pw_path(pw_id)
            print('Loading over/under expressed data {:3}/{} path={}'.format(ind+1, num_pw, path), end='\r')
            all_pw_map[pw_id] = nx.read_gpickle(path)
        print()
        return all_pw_map

    def save_pathways():
        for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
            path = get_pw_path(pw_id)
            print('Saving over/under expressed data {:3}/{} path={}'.format(ind+1, num_pw, path), end='\r')
            nx.write_gpickle(pw, path)
        print()
        return all_pw_map

    # check if we already stored all over/under expression pathway data if so restore them
    if np.all([os.path.exists(get_pw_path(pw_id)) for pw_id in all_pw_map]):
        return restore_pathways()

    num_pat = pat_ids.shape[0]
    # if there are missing ones calculate all of them
    print('Over and under expressed patient pathway labeling')
    for ind, pid in enumerate(pat_ids):
        print('Checking patient for over-expressed  {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
        gene_ind = (GE[..., pat_ids == pid] == 1).flatten() # over expressed genes
        genes = uni_ids[gene_ind] # get uniprot gene ids from indices
        label_mapper.mark_label_on_pathways('oe', pid, all_pw_map, genes, label)

        print('Checking patient for under-expressed {:4}/{} pid={}'.format(ind + 1, num_pat, pid))
        gene_ind = (GE[..., pat_ids == pid] == -1).flatten() # under expressed genes
        genes = uni_ids[gene_ind] # get uniprot gene ids from indices
        label_mapper.mark_label_on_pathways('ue', pid, all_pw_map, genes, label)

    return save_pathways()

label_patient_genes(all_pw_map, pat_ids, GE)

# experiment variables
smoothing_alpha = 0
a_smspk = smspk.smspk()

num_pat = pat_ids.shape[0]
num_pw = len(all_pw_map)
# calculate kernel matrices for over expressed genes
over_exp_kms = np.zeros((num_pw, num_pat, num_pat))
for ind, (pw_id, pw) in enumerate(all_pw_map.items()): # for each pathway
    over_exp_kms[ind] = a_smspk.kernel(pat_ids, pw, smoothing_alpha, label_key='label-oe', normalization=True)

# calculate kernel matrices for under expressed genes
under_exp_kms = np.zeros((num_pw, num_pat, num_pat))
for ind, (pw_id, pw) in enumerate(all_pw_map.items()): # for each pathway
    under_exp_kms[ind] = a_smspk.kernel(pat_ids, pw, smoothing_alpha,label_key='label-ue', normalization=True)


def plot_hm(data):
    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.show()

# pdb.set_trace()
print('End!')
