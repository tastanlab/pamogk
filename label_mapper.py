import numpy as np
import networkx as nx
import pdb

def mark_label_on_pathways(pid, pat_map, pw_map, gene_id_list, label=1):
	'''Marks given genes to the pathways

	Parameters
	----------
	pid: int
		patient id
	pw_map: map of networkx graphs of pathways
		patient label mapping
	pat_map: dict
		for a single patient; pathway to gene mapping
	gene_id_list: list of list of string
		uniprot gene id list of genes
	label: int
		the label which will be assigned to found genes in pathways - default value is 1
	'''
	gene_ids = [uid for a in gene_id_list for uid in a]
	for pw in pw_map.values(): # for each pathway
		for n in pw.nodes():
			nd = pw.nodes[n]
			if np.any([g in nd['uniprot_ids'] for g in gene_ids]):
				nd['label'][pid] = label
