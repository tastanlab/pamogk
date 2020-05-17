import argparse
import csv

from gensim.models.keyedvectors import KeyedVectors

from .. import config
from ..gene_mapper import uniprot_mapper as um
from ..pathway_reader import cx_pathway_reader

'''
Calling get_n2v_representations() returns vectors for each patient.
'''

KIRC_PATH = config.DATA_DIR / 'kirc_data' / 'kirc_somatic_mutation_data.csv'
UNIPROT_ENTREZ_MAP_FPATH = config.ROOT_DIR / 'gene_mapper' / 'uniprot-entrez-map.tsv'


def patient_entrez_to_uniprot():
    list_of_gene_patient = []
    with open(KIRC_PATH, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(csv_reader)
        for row in csv_reader:
            if int(row[1]) != 0:
                list_of_gene_patient.append(list(row))

    u_map = []
    e_map = []
    with open(UNIPROT_ENTREZ_MAP_FPATH, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        next(csv_reader)  # skip header
        for row in csv_reader:
            u_map.append(row[0])
            e_map.append(row[1])

    uni_to_entrez, entrez_to_uni = um.json_to_dict()

    patient_uniprot_list = []
    for row in list_of_gene_patient:
        uni_prot = []
        patient = row[2]

        try:
            uni_prot.append(entrez_to_uni[row[1]])
        except:
            if row[2] in e_map:
                uni_prot.append(u_map[e_map.index(row[2])])
            else:
                print("none")

        if len(uni_prot) != 0:
            patient_uniprot_list.append([patient, uni_prot])
    return patient_uniprot_list


def find_pathways_for_all_patients():
    patient_uniprot_list = patient_entrez_to_uniprot()
    count = 0
    pat_pws = []
    pathways = cx_pathway_reader.get_pathway_map()
    for patient in patient_uniprot_list:
        pathways_part = []
        count += 1
        if len(patient[1]) > 0:
            for geneList in patient[1]:
                x = cx_pathway_reader.get_pathways_with_genes(pathways, geneList)
                pathways_part.append(x)
        if len(pathways_part) != 0:
            pat_pws.append([patient[0], pathways_part])

    path = config.DATA_DIR / 'KircPathways.txt'
    with open(path, 'w') as f:
        for pat_pw in pat_pws:
            print(','.join(pat_pw[1]), file=f)


def get_patient_pathways_from_file():
    patient_pathway_list = []
    path = config.DATA_DIR / 'KircPathways.txt'
    with open(path, 'r') as f:
        for line in f:
            line = line[:-1]
            tokens = line.split(",")
            patient_pathway_list.append([tokens[0], tokens[1:]])
    return patient_pathway_list


def get_patient_uniprots_from_file():
    patient_uniprot_list = []
    path = config.DATA_DIR / 'KircUniprots.txt'
    with open(path, 'r') as f:
        for line in f:
            line = line[:-1]
            tokens = line.split(',')
            patient_uniprot_list.append([tokens[0], tokens[1:]])
    return patient_uniprot_list


def get_n2v_representations():
    parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
    parser.add_argument('--pathways', metavar='pathway-id', type=str, nargs='+', help='pathway ID list',
                        default=['hsa04151'])
    parser.add_argument('--debug', action='store_true', dest='debug', help='Enable Debug Mode')
    parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
    parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
    parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float,
                        help='Node2Vec feature space size', default=128)
    parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
    parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed',
                        default=False)
    parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments',
                        default=1000)
    parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float,
                        help='Surviving patient percentage in range [0, 1]', default=0.9)
    parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]',
                        default=0.4)

    args = parser.parse_args()
    patient_uniprot_list = get_patient_uniprots_from_file()
    patient_pathway_list = get_patient_pathways_from_file()
    count = 0
    patient_vector_list = []
    for patient_pathway in patient_pathway_list:
        vector_list = []
        for pathway in patient_pathway[1]:
            nx_g = cx_pathway_reader.read_single_pathway(pathway)
            # gene_vec_map = node2vec_processor.process(pathway, nx_g, args)
            filename = f'{pathway}-p={args.p:0.2f}-q={args.q:0.2f}-dir={args.is_directed}-run={args.rid}-word2vec.csv'
            gene_vec_map = {}
            filepath = config.DATA_DIR / 'node2vec' / filename
            key_vectors = KeyedVectors.load_word2vec_format(filepath, binary=False)
            for (eid, gene_vec) in zip(key_vectors.index2entity, key_vectors.vectors):
                gene_vec_map[int(eid)] = gene_vec
            uniprot_list = patient_uniprot_list[count][1]
            for gene in list(nx_g.nodes(data=True)):
                for alias in gene[1]['alias']:
                    if len(alias.split(':')) > 1:
                        alias = alias.split(':')[1]
                    if alias in uniprot_list:
                        vector_list.append(gene_vec_map[gene[0]])
        patient_vector_list.append([patient_pathway[0], vector_list])
        count += 1

    path = config.DATA_DIR / 'KircVectors.txt'
    with open(path, 'w') as f:
        for patient in patient_vector_list:
            line = patient[0]
            for vector in patient[1]:
                line += ";" + ','.join(map(str, vector))
            f.write(line + "\n")
    return patient_vector_list

# get_n2v_representations()
# find_pathways_for_all_patients()
