#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
import scipy
import pdb

class smspk:

	@staticmethod
	def kernel(data, alpha, epsilon=10**-6, normalization=False):
		# data: a list of networkx graphs
		# alpha: the smoothing parameter
		# epsilon: smoothing converges if the change is lower than epsilon -- default value is 10^-6
		# normalization: normalize the kernel matrix such that the diagonal is 1 -- default value is False

		# extract labels of nodes of graphs -- ASSUMPTION: all graphs have the same nodes
		nodes = list(nx.get_node_attributes(data[0], 'label').keys())
		mutations = np.empty([len(data), len(nodes)])
		for i in range(len(data)):
			tmp = nx.get_node_attributes(data[i], 'label')
			for j in range(len(nodes)):
				mutations[i,j] = tmp[nodes[j]]

		# extract the adjacency matrix on the order of nodes we have
		adj_mat = nx.to_numpy_array(data[0], nodelist=nodes)

		# smooth the mutations through the pathway
		mutations = smspk.smooth(mutations, adj_mat, alpha, epsilon)

		# find shortest paths between all pairs of nodes which are genes
		all_sp = nx.all_pairs_shortest_path(data[0])
		node_types = nx.get_node_attributes(data[0], 'type')

		km = np.zeros((len(data), len(data)))
		# print km
		skip = 1 # it is used to indicate how many nodes from the beginning will be skipped in the shortest path list -- we do not want to process the same shortest paths again
		for a_sp in all_sp:
			if node_types[a_sp[0]] == 'Protein': # if the source is gene/protein
				tmp_sp_of_nodes = list(a_sp[1].keys())
				for i in range(skip, len(tmp_sp_of_nodes)):
					if node_types[a_sp[1][tmp_sp_of_nodes[i]][-1]] == 'Protein': # if the destination is gene/protein
						# print("Shortest path: {}".format(a_sp[1][tmp_sp_of_nodes[i]]))
						ind = np.isin(nodes, a_sp[1][tmp_sp_of_nodes[i]])
						tmp_md = mutations[:][:,ind]
						tmp_km = np.matmul(tmp_md, np.transpose(tmp_md)) # calculate similarities of patients based on the current pathway
						km += tmp_km # update the main kernel matrix
			skip += 1

		# normalize the kernel matrix if normalization is true
		if normalization == True:
			km = smspk.normalize_kernel_matrix(km)

		return km


	@staticmethod
	def smooth(md, adj_m, alpha, epsilon=10**-6):
		# md: a numpy array of genes of patients indicating which one is mutated or not
		# adj_m: the adjacency matrix of the pathway
		# alpha: the smoothing parameter
		# epsilon: smoothing converges if the change is lower than epsilon -- default value is 10^-6
		norm_adj_mat = adj_m @ np.diag(1.0 / np.sum(adj_m, axis=0))

		s_md = md
		pre_s_md = md + epsilon + 1

		while np.linalg.norm(s_md - pre_s_md) > epsilon:
			pre_s_md = s_md
			s_md = ((alpha * pre_s_md) @ norm_adj_mat) + (1 - alpha) * md

		return s_md

	@staticmethod
	def normalize_kernel_matrix(km):
		D = np.diag(1 / np.sqrt(np.diag(km)))
		norm_km = D @ km @ D # K_ij / sqrt(K_ii * K_jj)
		return np.nan_to_num(norm_km) # replace NaN with 0

def main():
	# Create a networkx graph object
	mg = nx.Graph()

	# Add edges to to the graph object
	# Each tuple represents an edge between two nodes
	mg.add_edges_from([
							('A','B'),
							('A','C'),
							('C','D'),
							('A','E'),
							('C','E'),
							('D','B'),
							('B','C'),
							('C','X')])

	mg2 = nx.Graph()
	mg2.add_edges_from([
							('A','B'),
							('A','C'),
							('C','D'),
							('A','E'),
							('C','E'),
							('D','B'),
							('B','C'),
							('C','X')])

	nx.set_node_attributes(mg, {'X':0, 'A':1, 'B':0, 'C':0, 'D':1, 'E':0}, 'label')
	#nx.set_node_attributes(mg, {'X':'Protein', 'A':'Calcium', 'B':'Protein', 'C':'Protein', 'D':'Calcium', 'E':'Protein'}, 'type')
	nx.set_node_attributes(mg, {'X':'Protein', 'A':'Protein', 'B':'Protein', 'C':'Protein', 'D':'Protein', 'E':'Protein'}, 'type')

	nx.set_node_attributes(mg2, {'X':1, 'A':1, 'B':0, 'C':0, 'D':1, 'E':0}, 'label')
	#nx.set_node_attributes(mg2, {'X':'Protein', 'A':'Calcium', 'B':'Protein', 'C':'Protein', 'D':'Calcium', 'E':'Protein'}, 'type')
	nx.set_node_attributes(mg2, {'X':'Protein', 'A':'Protein', 'B':'Protein', 'C':'Protein', 'D':'Protein', 'E':'Protein'}, 'type')

	# smoothing parameter for the kernel
	alpha = 0.1

	# calculate the kernel using smspk
	km = smspk.kernel([mg, mg2], alpha)

	# display the resulting kernel matrix
	print(("Kernel matrix calculated by smspk with alpha {}:".format(alpha)))
	print(km)



if __name__ == '__main__':
	main()


#atadam
