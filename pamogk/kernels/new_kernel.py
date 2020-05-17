import argparse

import pandas

from pamogk import config
from ..gene_mapper import uniprot_mapper
from ..kernels import node2vec_h_i_k as n2v
from ..lib.sutils import *
from ..pathway_reader import cx_pathway_reader as cx_pw

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--patient-data', '-f', metavar='file-path', dest='patient_data', type=Path, help='pathway ID list',
                    default=config.DATA_DIR / 'kirc_data/kirc_somatic_mutation_data.csv')

args = parser.parse_args()
print_args(args)

p, q = 1, 1
num_walks, walk_length = 10, 20
output = "xd.emb"


@timeit
def read_data():
    # Real Data #
    # process RNA-seq expression data
    patients = {}
    with open(args.patient_data) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pat_id = row['Patient ID']
            ent_id = row['Entrez Gene ID']
            if pat_id not in patients:
                patients[pat_id] = {ent_id}
            else:
                patients[pat_id].add(ent_id)

    return patients


@timeit
def preprocess_patient_data(patients):
    # get the dictionary of gene id mappers
    uni2ent, ent2uni = uniprot_mapper.json_to_dict()

    res = []
    num_empty = 0
    for pat_id, ent_ids in patients.items():
        # uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
        uni_ids = [uid for eid in ent_ids if eid in ent2uni for uid in ent2uni[eid]]
        # if there are any matches map them
        if len(uni_ids) > 0:
            res.append({
                'pat_id': pat_id,
                'mutated_nodes': uni_ids,
            })
        else:
            num_empty += 1

    log('removed patients:', num_empty)

    return res


def eliminate_with_conf(walk, conf, nodes):
    nodeFreq = np.zeros(len(nodes))
    len_walk = len(walk)
    threshold = conf * len_walk
    for w in walk:
        for nodeId in w:
            nodeFreq[nodeId] += 1

    return [idx for idx, freq in enumerate(nodeFreq) if freq >= threshold]


def process_walks(walks, conf, nodes):
    for node in nodes:
        walks[node] = eliminate_with_conf(walks[node], conf, nodes)


def get_neighbors_in_pathway(pw_graph, conf):
    nx_G = pw_graph

    id_gene_map = {}
    for nodeId in pw_graph._node.keys():
        id_gene_map[nodeId] = pw_graph._node[nodeId]["uniprotids"]

    directed = False
    G = n2v.Graph(nx_G, directed, p, q)
    G.preprocess_transition_probs()
    walks = G.simulate_walks(num_walks, walk_length)
    process_walks(walks, conf, walks.keys())
    # np.savetxt(output, walks)
    return walks, id_gene_map


@timeit
def read_pathways():
    # get all pathways
    return cx_pw.read_pathways()


@timeit
def get_neighbors_for_all_pathways(all_pw_map, conf):
    res = {}
    map = {}
    for ind, (pw_id, pw) in enumerate(all_pw_map.items()):
        neighbors, mapper = get_neighbors_in_pathway(pw, conf)
        res[pw_id] = neighbors
        map[pw_id] = mapper
        break
    return res, map


def calc_patientwise_score(neighbors, patient1, patient2, mapper):
    # from patient1 to patient2
    # find genes in graph and find unique set
    patient1_nodes = []
    for node in neighbors.keys():
        if len(set(mapper[node]).intersection(patient1["mutated_nodes"])) > 0:
            patient1_nodes.append(node)

    hit_count = 0
    for p1_node in patient1_nodes:
        for neigh in neighbors[p1_node]:
            if len(set(mapper[neigh]).intersection(patient2["mutated_nodes"])) > 0:
                hit_count += 1
                break
    if len(patient1_nodes) == 0:
        p1_rate = 0
    else:
        p1_rate = float(hit_count) / len(patient1_nodes)

    # from patient2 to patient1
    patient2_nodes = []
    for node in neighbors.keys():
        if len(set(mapper[node]).intersection(patient2["mutated_nodes"])) > 0:
            patient2_nodes.append(node)

    hit_count = 0
    for p2_node in patient2_nodes:
        for neigh in neighbors[p2_node]:
            if len(set(mapper[neigh]).intersection(patient1["mutated_nodes"])) > 0:
                hit_count += 1
                break

    if len(patient2_nodes) == 0:
        p2_rate = 0
    else:
        p2_rate = float(hit_count) / len(patient2_nodes)

    avg = (p1_rate + p2_rate) / 2

    return avg


@timeit
def calc_similarity_from_pathway(neighbors, patients, id_mapper):
    len_p = len(patients)
    similarityMatrix = np.zeros((len_p, len_p))
    for i in range(len_p):
        for j in range(i, len_p):
            score = calc_patientwise_score(neighbors, patients[i], patients[j], id_mapper)
            if i == j and score != 1:
                similarityMatrix[i, j] = 1
            else:
                similarityMatrix[i, j] = score
                similarityMatrix[j, i] = score

    return similarityMatrix


@timeit
def calc_kernel_from_similarity(similarityMatrix):
    distanceMatrix = 1 - similarityMatrix
    sigmasqList = [0.2, 0.5, 1, 5]
    sigmasq = sigmasqList[1]
    # it was scipt.exp but they should be doing the same thing
    return np.exp(-np.square(distanceMatrix) / (2 * sigmasq))


@timeit
def calc_kernel_from_pathways(neighbor_mappings, patients, id_mapper):
    kernel_res = None
    flag = 0
    len_p = len(patients)
    for pathway_id in neighbor_mappings.keys():
        similarityMatrix = calc_similarity_from_pathway(neighbor_mappings[pathway_id], patients, id_mapper[pathway_id])
        one_kernel = calc_kernel_from_similarity(similarityMatrix)
        if flag == 0:
            kernel_res = one_kernel
            kernel_res = np.reshape(kernel_res, (-1, len_p, len_p))
            flag = 1
        else:
            one_kernel = np.reshape(one_kernel, (-1, len_p, len_p))
            kernel_res = np.concatenate((kernel_res, one_kernel), axis=0)

    return kernel_res


def load_kernel(file_loc):
    data = pandas.read_csv(file_loc, delimiter=" ", dtype=np.float64)
    return np.array(data).reshape((165, 417, 417))


def isPSD(A, tol=1e-8):
    E, V = np.linalg.eigh(A)
    return np.all(E >= -tol)
