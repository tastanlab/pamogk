import sys
import numpy as np
from gene_mapper import uniprot_mapper as um
from pathway_reader import cx_pathway_reader as cx_pw
import rnaseq_processor as rp
from operator import itemgetter
import operator
import networkx as nx
import pdb

def map_label_on_pathways(pw_list, gene_label_list, label=1):
	'''
	Inputs
		pw_list: a list of networkx instances of pathways
		gene_id: uniprot id of the gene
		label: the label which will be assigned to found genes in pathways - default value is 1

	Outputs:
		pw_list: a list of modified (i.e. label of gene is mapped to corresponding nodes) networkx instances of pathways
	'''
	# print(type(gene_label_list))
	# print(gene_label_list)
	gene_label_dict = {}
	whole_genes = [item for sublist in gene_label_list for item in sublist]

	# pdb.set_trace()

	for i in pw_list.keys(): # for each pathway
		alias = nx.get_node_attributes(pw_list[i], 'alias')
		tmp_labels = dict(zip(alias.keys(), np.zeros(len(alias), dtype=int)))
		for j in alias.keys():
			for t in range(len(alias[j])):
				if alias[j][t].split(':')[1] in whole_genes:
					tmp_labels[j] = label
					break
		nx.set_node_attributes(pw_list[i], tmp_labels, 'label')

	return pw_list

	# for key in pw_dict:
	# 	nodes = pw_dict[i].nodes()
	# 	tmp_dict = dict(zip(nodes, np.zeros((len(nodes),), dtype=int)))
	# 	tmp_dict.update(gene_label_dict)
	# 	nx.set_node_attributes(pw_dict[i], tmp_dict, 'label')
    #
	# return pw_dict












#atadam
