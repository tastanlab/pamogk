import networkx as nx
import numpy as np
import scipy

class smspk:

    @staticmethod
    def kernel(data, alpha, epsilon=10**-6):
        # data: a list of networkx graphs
        # alpha: the smoothing parameter
        # epsilon: smoothing converges if the change is lower than epsilon -- default value is 10^-6

        # extract labels of nodes of graphs -- ASSUMPTION: all graphs have the same nodes
        nodes = nx.get_node_attributes(data[0], 'label').keys()
        mutations = np.empty([len(data), len(nodes)])
        for i in xrange(len(data)):
            tmp = nx.get_node_attributes(data[i], 'label')
            for j in xrange(len(nodes)):
                mutations[i,j] = tmp[j]

        # extract the adjacency matrix on the order of nodes we have
        adj_mat = nx.to_numpy_array(data[0], nodelist=nodes)

        # smooth the mutations through the pathway
        mutations = smspk.smooth(mutations, adj_mat, alpha, epsilon)

        # find shortest paths between all pairs of nodes
        all_sp = nx.all_pairs_shortest_path(data[0])

        km = np.empty([len(data), len(data)])
        skip = 1 # it is used to indicate how many nodes from the beginning will be skipped in the shortest path list -- we do not want to process the same shortest paths again
        for a_sp in all_sp:
            tmp_sp_of_nodes = a_sp[1].keys()
            for i in xrange(skip, len(tmp_sp_of_nodes)):
                ind = np.isin(nodes, a_sp[1][tmp_sp_of_nodes[i]])
                tmp_md = mutations[:][:,ind]
                tmp_km = np.matmul(tmp_md, np.transpose(tmp_md)) # calculate similarities of patients based on the current pathway
                km += tmp_km # update the main kernel matrix
            skip += 1

        return km


    @staticmethod
    def smooth(md, adj_m, alpha, epsilon=10**-6):
            # md: a numpy array of genes of patients indicating which one is mutated or not
            # adj_m: the adjacency matrix of the pathway
            # alpha: the smoothing parameter
            # epsilon: smoothing converges if the change is lower than epsilon -- default value is 10^-6

            norm_adj_mat = np.matmul(adj_m, np.diag(1.0 / np.sum(adj_m, axis=0)))

            s_md = md
            pre_s_md = md + epsilon + 1

            while np.linalg.norm(s_md - pre_s_md) > epsilon:
                pre_s_md = s_md
                s_md = np.matmul(alpha * pre_s_md, norm_adj_mat) + (1 - alpha) * md

            return s_md





#atadam
