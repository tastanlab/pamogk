import numpy as np
import networkx as nx
import math
import random
import sys

hsa_04151 = {
    'cell_cycle_nodes': [24, 35, 105, 62, 26, 25, 61],
    'cell_survival_nodes': [27, 28, 32, 98, 96, 103, 110, 97, 99, 33],
    'start_nodes': [43, 86, 78, 91, 104, 84],
}


def generate_patients(G, num_pat, surv_dist, mut_dist=0.2, psm=0.9, has_cycle=False, has_survival=True):
    print('Generating %d patients with %2d%% survival rate:' % (num_pat, surv_dist * 100))
    patients = []
    num_surv = math.ceil(surv_dist * num_pat) # surviving patient count
    for i in range(num_pat):
        sick = i >= num_surv
        patient = { 'pid': i, 'mutated_nodes': [] }
        patient['sick'] = sick
        _psm = psm if sick else 1 - psm
        patient['mutated_nodes'] = calc_mutated_nodes(G, mut_dist, _psm, has_cycle, has_survival)
        patients.append(patient)
    random.shuffle(patients)
    for i in range(num_pat):
        patients[i]['pid'] = i
    return patients
'''
This is for hsa04151, for others we don't have same structure
'''
def calc_mutated_nodes(G, mut_dist=0.2, psm=0.9, has_cycle=False, has_survival=True):
    dest_nodes = []
    if has_cycle: dest_nodes += hsa_04151['cell_cycle_nodes']
    if has_survival: dest_nodes += hsa_04151['cell_survival_nodes']
    surv_nodes = set()
    for st in hsa_04151['start_nodes']:
        paths = nx.single_source_shortest_path(G, st)
        # print(paths)
        for dt in dest_nodes:
            if dt in paths:
                p = paths[dt]
                print('Found path (%3d, %3d):' % (st, dt), p)
                surv_nodes |= set(p)
    print('Survival nodes:', surv_nodes)
    print('Probability of Survival Node is set to:', psm)
    print('Mutation Distribution is set to:', mut_dist, 'of all genes')

    num_mut_nodes = math.ceil(len(G.nodes()) * mut_dist)
    num_srv_mut_nodes = math.ceil(num_mut_nodes * psm)
    print('Nodes to mutate:', num_mut_nodes)
    print('Survival Nodes to mutate:', num_srv_mut_nodes)

    seed = random.randrange(sys.maxsize)
    rand = random.Random(seed)
    print('Setting random seed to:', seed)

    # mutated nodes to return
    mutated_nodes = set()
    print('Assigning mutated survival nodes')
    mutated_nodes |= mutate_nodes(surv_nodes, num_srv_mut_nodes, rand)
    print('Assigning mutated non-survival nodes')
    # get non survival nodes
    non_surv_nodes = [eid for eid in G.nodes() if eid not in surv_nodes]
    mutated_nodes |= mutate_nodes(non_surv_nodes, num_mut_nodes - num_srv_mut_nodes, rand)

    print('All mutated nodes:', mutated_nodes)
    return mutated_nodes

def mutate_nodes(node_list, node_mut_count, rand):
    mutated_nodes = set()
    # copy to not mutate original
    nlist = list(node_list)
    for i in range(node_mut_count):
        len_nodes = len(nlist)
        if len_nodes == 0:
            print('Warning! Not enough nodes for mutation! Could not assign:', node_mut_count - i - 1)
            break
        mut_i = rand.randint(0, len(nlist) - 1)
        mut_n = nlist.pop(mut_i)
        mutated_nodes.add(mut_n)
    print('Mutated nodes:', mutated_nodes)
    return mutated_nodes

