import numpy as np
import networkx as nx
import pandas as pd
import copy
import pdb
from operator import itemgetter
import operator
import matplotlib.pyplot as plt
import smspk
import label_mapper
import rnaseq_processor as rp
from pathway_reader import cx_pathway_reader as cx_pw
from gene_mapper import uniprot_mapper as um

### Real Data ###
# process RNA-seq expression data
[gene_exp, genes_map] = rp.process("/home/aburak/Projects/smSPK/research_codes/data/kidney_files/unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt")
# pdb.set_trace()

# get the dictionary of gene id mappers
[uni2ent, ent2uni] = um.jsonToDict()

# convert entrez gene id to uniprot id
tmp = gene_exp.index.tolist()
uni_of_rna = [None] * len(tmp)
for i in range(len(tmp)):
	if tmp[i] in ent2uni:
		uni_of_rna[i] = ent2uni[tmp[i]]
# pdb.set_trace()
none_ind = np.where(np.array(uni_of_rna) == None)[0]
print("Number of genes whose uniprot id is not found: ", len(none_ind))

# prune genes whose uniprot id is not found
gene_exp = gene_exp.drop(gene_exp.index[none_ind])
pruned_uni_of_rna = itemgetter(*(set(range(len(uni_of_rna))).difference(set(none_ind))))(uni_of_rna)

gene_exp["uniprot_id"] = pruned_uni_of_rna
print(len(pruned_uni_of_rna))

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
# pathways = {"pw1": mg, "pw2": mg2}

# pdb.set_trace()
###############################################################################
# experiment variables
smoothing_alpha = 0
a_smspk = smspk.smspk()
###############################################################################
# over-expressed genes
print("Over-expressed genes will be used to calculate kernel matrices...")
over_exp_pathways_of_patients = dict(zip(gene_exp.columns[0:-1], [copy.deepcopy(pathways) for i in range(len(gene_exp.columns)-1)]))
for i in over_exp_pathways_of_patients.keys():
	print(i)
	tmp = gene_exp.index[gene_exp[i] == 1]
	over_exp_pathways_of_patients[i] = label_mapper.map_label_on_pathways(over_exp_pathways_of_patients[i], gene_exp.loc[tmp,"uniprot_id"].tolist())

# calculate kernel matrices
a_patient_key = list(over_exp_pathways_of_patients.keys())[0] # to get the size of pathways
over_exp_kms = np.zeros((len(over_exp_pathways_of_patients[a_patient_key]), len(over_exp_pathways_of_patients), len(over_exp_pathways_of_patients)))
ind = 0
for pw_key in over_exp_pathways_of_patients[a_patient_key].keys(): # for each pathway
	tmp = [over_exp_pathways_of_patients[ptnt_key][pw_key] for ptnt_key in over_exp_pathways_of_patients.keys()] # list of the same pathway of all patients
	over_exp_kms[ind] = a_smspk.kernel(tmp, smoothing_alpha, normalization=True)
	ind += 1

del over_exp_pathways_of_patients

###############################################################################
# under-expressed genes
print("Under-expressed genes will be used to calculate kernel matrices...")
under_exp_pathways_of_patients = dict(zip(gene_exp.columns[0:-1], [copy.deepcopy(pathways) for i in range(len(gene_exp.columns)-1)]))
for i in under_exp_pathways_of_patients.keys():
	print(i)
	tmp = gene_exp.index[gene_exp[i] == -1]
	under_exp_pathways_of_patients[i] = label_mapper.map_label_on_pathways(under_exp_pathways_of_patients[i], gene_exp.loc[tmp,"uniprot_id"].tolist())

# calculate kernel matrices
a_patient_key = list(under_exp_pathways_of_patients.keys())[0] # to get the size of pathways
under_exp_kms = np.zeros((len(under_exp_pathways_of_patients[a_patient_key]), len(under_exp_pathways_of_patients), len(under_exp_pathways_of_patients)))
ind = 0
smoothing_alpha = 0
for pw_key in under_exp_pathways_of_patients[a_patient_key].keys(): # for each pathway
	tmp = [under_exp_pathways_of_patients[ptnt_key][pw_key] for ptnt_key in under_exp_pathways_of_patients.keys()] # list of the same pathway of all patients
	under_exp_kms[ind] = a_smspk.kernel(tmp, smoothing_alpha, normalization=True)
	ind += 1

###############################################################################


def plot_hm(data):
	plt.imshow(data, cmap='hot', interpolation='nearest')
	plt.show()

# pdb.set_trace()
print("End!")









# atadam
