import numpy as np
import networkx as nx

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
	whole_genes = [item for sublist in gene_label_list for item in sublist]

	for pw in pw_list.values(): # for each pathway
		alias = nx.get_node_attributes(pw, 'alias')
		tmp_labels = dict(zip(alias.keys(), np.zeros(len(alias), dtype=int)))
		for j in alias.keys():
			for t in range(len(alias[j])):
				if alias[j][t].split(':')[1] in whole_genes:
					tmp_labels[j] = label
					break
		nx.set_node_attributes(pw, tmp_labels, 'label')
