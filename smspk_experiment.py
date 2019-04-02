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
pathways = cx_pw.read_pathways()

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
# pathways = {'pw1': mg, 'pw2': mg2}

class timeit(object):
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        print("Started:", self.f.__name__)
        t = time.time()
        res = self.f(*args, **kwargs)
        print("Finished:", self.f.__name__, " elapsed:", int(time.time() - t))
        return res

@timeit
def clone_pathway_map(pathways, cols):
	return dict((c, dict((k, v.copy()) for k, v in pathways.items())) for c in cols)

# pdb.set_trace()
###############################################################################
# experiment variables
smoothing_alpha = 0
a_smspk = smspk.smspk()
###############################################################################
# over-expressed genes
print('Over-expressed genes will be used to calculate kernel matrices...')
OEPP = clone_pathway_map(pathways, pat_ids) # over expressed pathways of patients
for pid, pat_pw_list in OEPP.items():
	print('Checking patient:', pid)
	pdb.set_trace()
	label_mapper.map_label_on_pathways(pat_pw_list, GE[..., pat_ids == pid].flatten())

# calculate kernel matrices
a_patient_key = list(OEPP.keys())[0] # to get the size of pathways
over_exp_kms = np.zeros((len(OEPP[a_patient_key]), len(OEPP), len(OEPP)))
ind = 0
for pw_key in OEPP[a_patient_key].keys(): # for each pathway
	tmp = [OEPP[ptnt_key][pw_key] for ptnt_key in OEPP.keys()] # list of the same pathway of all patients
	over_exp_kms[ind] = a_smspk.kernel(tmp, smoothing_alpha, normalization=True)
	ind += 1

del OEPP

###############################################################################
# under-expressed genes
print('Under-expressed genes will be used to calculate kernel matrices...')
UEPP = clone_pathway_map(pathways, pat_ids) # under expressed pathways of patients
for pid, pat_pw_list in OEPP.items():
	print('Checking patient under-epxressed:', pid)
	tmp = gene_exp.index[gene_exp[pid] == -1]
	label_mapper.map_label_on_pathways(pat_pw_list, GE[..., pat_ids == pid].flatten())

# calculate kernel matrices
a_patient_key = list(UEPP.keys())[0] # to get the size of pathways
under_exp_kms = np.zeros((len(UEPP[a_patient_key]), len(UEPP), len(UEPP)))
ind = 0
smoothing_alpha = 0
for pw_key in UEPP[a_patient_key].keys(): # for each pathway
	tmp = [UEPP[ptnt_key][pw_key] for ptnt_key in UEPP.keys()] # list of the same pathway of all patients
	under_exp_kms[ind] = a_smspk.kernel(tmp, smoothing_alpha, normalization=True)
	ind += 1

###############################################################################


def plot_hm(data):
	plt.imshow(data, cmap='hot', interpolation='nearest')
	plt.show()

# pdb.set_trace()
print('End!')









# atadam
