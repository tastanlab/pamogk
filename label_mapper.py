import numpy as np
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
	gene_label_dict = {}

	for pw in pw_list.values(): # for each pathway
		pdb.set_trace()
		alias_list = [a.split(':')[1] for a in nx.get_node_attributes(pw, 'alias')]
		nx.set_node_attributes(pw, [label if a in gene_label_list else 0 for a in alias_list], 'label')
