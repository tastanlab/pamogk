import networkx as nx
from pamogk.lib.sutils import *


# options are raise, warn, print (default)
# np.seterr(all='raise')

def kernel(pat_ids, pathway, label_key, alpha=0.5, epsilon=1e-6, normalization=False):
    """
    Parameters
    ----------
    pat_ids:
        list of patient ids
    pathway:
        pathway networkx graph
    alpha: float
        the smoothing parameter
    label_key: str
        label attribute key showing patient to label mapping. Should be filled in experiment
    epsilon: {1e-6} float
        smoothing converges if the change is lower than epsilon
    normalization: {False} bool
        normalize the kernel matrix such that the diagonal is 1
    """

    num_pat = pat_ids.shape[0]
    pat_ind = {}
    for ind, pid in enumerate(pat_ids): pat_ind[pid] = ind
    # extract labels of nodes of graphs
    mutations = np.zeros([num_pat, len(pathway.nodes)])
    for nid in pathway.nodes:
        nd = pathway.nodes[nid]
        for pid, lb in nd[label_key].items():
            if pid in pat_ind.keys():
                mutations[pat_ind[pid], nid] = lb

    # extract the adjacency matrix on the order of nodes we have
    adj_mat = nx.to_numpy_array(pathway, nodelist=pathway.nodes)
    ordered_graph = nx.OrderedGraph()
    ordered_graph.add_nodes_from(pathway.nodes())
    ordered_graph.add_edges_from(sorted(list(pathway.edges())))

    # smooth the mutations through the pathway
    mutations = smooth(mutations, adj_mat, alpha, epsilon)
    # get all pairs shortest paths
    all_pairs_sp = nx.all_pairs_shortest_path(ordered_graph)

    km = np.zeros((num_pat, num_pat))

    checked = []
    for src, dsp in all_pairs_sp:  # iterate all pairs shortest paths
        # add source node to checked nodes so we won't check it again in destinations
        checked.append(src)
        # skip if the source is not gene/protein
        if pathway.nodes[src]['type'] != 'Protein': continue
        # otherwise
        for dst, sp in dsp.items():
            # if destination already checked skip
            if dst in checked: continue
            # if the destination is not gene/protein skip
            if pathway.nodes[sp[-1]]['type'] != 'Protein': continue
            ind = np.isin(pathway.nodes, sp)
            tmp_md = mutations[:, ind]
            # calculate similarities of patients based on the current pathway
            tmp_km = tmp_md @ np.transpose(tmp_md)
            km += tmp_km  # update the main kernel matrix

    # normalize the kernel matrix if normalization is true
    if normalization == True: km = normalize_kernel_matrix(km)

    return km


def smooth(md, adj_m, alpha=0.5, epsilon=10 ** -6):
    """
    md: numpy array
        a numpy array of genes of patients indicating which one is mutated or not
    adj_m: numpy array
        the adjacency matrix of the pathway
    alpha: {0.5} float
        the smoothing parameter in range of 0-1
    epsilon: {1e-6} float
        smoothing converges if the change is lower than epsilon
    """
    # since alpha will be together with norm_adj_mat all the time multiply here
    alpha_norm_adj_mat = alpha * adj_m / np.sum(adj_m, axis=0)

    s_md = md
    pre_s_md = md + epsilon + 1

    while np.linalg.norm(s_md - pre_s_md) > epsilon:
        pre_s_md = s_md
        # alpha_norm_adj_mat already inclues alpha multiplier
        s_md = (s_md @ alpha_norm_adj_mat) + (1 - alpha) * md

    return s_md


def normalize_kernel_matrix(km):
    kmD = np.array(np.diag(km))
    kmD[kmD == 0] = 1
    D = np.diag(1 / np.sqrt(kmD))
    norm_km = D @ km @ D  # K_ij / sqrt(K_ii * K_jj)
    return np.nan_to_num(norm_km)  # replace NaN with 0


def main():
    # Create a networkx graph object
    mg = nx.Graph()

    # Add edges to to the graph object
    # Each tuple represents an edge between two nodes
    mg.add_edges_from([
        ('A', 'B'),
        ('A', 'C'),
        ('C', 'D'),
        ('A', 'E'),
        ('C', 'E'),
        ('D', 'B'),
        ('B', 'C'),
        ('C', 'X')])

    mg2 = nx.Graph()
    mg2.add_edges_from([
        ('A', 'B'),
        ('A', 'C'),
        ('C', 'D'),
        ('A', 'E'),
        ('C', 'E'),
        ('D', 'B'),
        ('B', 'C'),
        ('C', 'X')])

    nx.set_node_attributes(mg, {'X': 0, 'A': 1, 'B': 0, 'C': 0, 'D': 1, 'E': 0}, 'label')
    # nx.set_node_attributes(mg, {'X':'Protein', 'A':'Calcium', 'B':'Protein', 'C':'Protein', 'D':'Calcium', 'E':'Protein'}, 'type')
    nx.set_node_attributes(mg, {'X': 'Protein', 'A': 'Protein', 'B': 'Protein', 'C': 'Protein', 'D': 'Protein',
                                'E': 'Protein'}, 'type')

    nx.set_node_attributes(mg2, {'X': 1, 'A': 1, 'B': 0, 'C': 0, 'D': 1, 'E': 0}, 'label')
    # nx.set_node_attributes(mg2, {'X':'Protein', 'A':'Calcium', 'B':'Protein', 'C':'Protein', 'D':'Calcium', 'E':'Protein'}, 'type')
    nx.set_node_attributes(mg2, {'X': 'Protein', 'A': 'Protein', 'B': 'Protein', 'C': 'Protein', 'D': 'Protein',
                                 'E': 'Protein'}, 'type')

    # smoothing parameter for the kernel
    alpha = 0.1

    # calculate the kernel using PAMOGK
    # NOTE: this might not be working after some changes to the kernel method
    km = kernel(np.array([0, 1]), mg, alpha=alpha, label_key='label')

    # display the resulting kernel matrix
    print('Kernel matrix calculated by PAMOGK with alpha', alpha)
    print(km)


if __name__ == '__main__':
    main()
