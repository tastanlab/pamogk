import sys
import config
sys.path.append(config.root_dir)
import os
import csv
from gene_mapper import uniprot_mapper as um
from pathway_reader import cx_pathway_reader
import argparse
from gensim.models.keyedvectors import KeyedVectors

'''
Calling get_n2v_representations() returns vectors for each patient.
'''

KIRC_PATH = os.path.join(config.data_dir, 'kirc_data', 'kirc_somatic_mutation_data.csv')
UNIPROT_ENTREZ_MAP_FPATH = os.path.join(config.root_dir, 'gene_mapper', 'uniprot-entrez-map.tsv')


def patient_entrez_to_uniprot():

    listOfGenePatient = []
    with open(KIRC_PATH, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(spamreader)
        for row in spamreader:
            if int(row[1])!=0:
                listOfGenePatient.append(list(row))

    UMap = []
    EMap = []
    with open(UNIPROT_ENTREZ_MAP_FPATH, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        next(spamreader)
        for row in spamreader:
            UMap.append(row[0])
            EMap.append(row[1])

    uniToEntrez, entrezToUni = um.json_to_dict()

    patient_uniprot_list = []
    for row in listOfGenePatient:
        uniProt=[]
        patient = row[2]

        try:
            uniProt.append(entrezToUni[row[1]])
        except:
            if row[2] in EMap:
                uniProt.append(UMap[EMap.index(row[2])])
            else:
                print("none")

        if len(uniProt)!=0:
            patient_uniprot_list.append([patient,uniProt])
    return patient_uniprot_list


def findPathwaysForAllPatients():
    patient_uniprot_list = patient_entrez_to_uniprot()
    count = 0
    pat_pws = []
    pathways = cx_pathway_reader.get_pathway_map()
    for patient in patient_uniprot_list:
        pathwaysPart = []
        count+=1
        if len(patient[1])>0:
            for geneList in patient[1]:
                x = cx_pathway_reader.get_pathways_with_gene(pathways, geneList)
                pathwaysPart.append(x)
        if len(pathwaysPart) != 0:
            pat_pws.append([patient[0], pathwaysPart])

    path = os.path.join(config.data_dir, 'KircPathways.txt')
    with open(path, 'w') as f:
        for pat_pw in pat_pws:
            print(','.join(pat_pw[1]), file=f)

def get_patient_pathways_from_file():
    patient_Pathway_List = []
    path = os.path.join(config.data_dir, 'KircPathways.txt')
    with open(path, 'r') as f:
        for line in f:
            line = line[:-1]
            tokens = line.split(",")
            patient_Pathway_List.append([tokens[0],tokens[1:]])
    return patient_Pathway_List


def get_patient_uniprots_from_file():
    patient_uniprot_list = []
    path = os.path.join(config.data_dir, 'KircUniprots.txt')
    with open(path, 'r') as f:
        for line in f:
            line = line[:-1]
            tokens = line.split(',')
            patient_uniprot_list.append([tokens[0],tokens[1:]])
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
    patient_Pathway_List = get_patient_pathways_from_file()
    count = 0
    patient_vector_list = []
    for patientPathway in patient_Pathway_List:
        vectorList = []
        for pathway in patientPathway[1]:
            nx_G = cx_pathway_reader.read_single_pathway(pathway)
            #gene_vec_map = node2vec_processor.process(pathway, nx_G, args)
            FNAME = '{}-p={:0.2f}-q={:0.2f}-dir={}-run={}-word2vec.csv'\
                .format(pathway, args.p, args.q, args.is_directed, args.rid)
            gene_vec_map = {}
            FPATH = os.path.join(config.data_dir, 'node2vec', FNAME)
            keyVectors = KeyedVectors.load_word2vec_format(FPATH, binary=False)
            for (eid, gene_vec) in zip(keyVectors.index2entity, keyVectors.vectors):
                gene_vec_map[int(eid)] = gene_vec
            uniprotList = patient_uniprot_list[count][1]
            for gene in list(nx_G.nodes(data=True)):
                for alias in gene[1]['alias']:
                    if len(alias.split(':')) > 1:
                        alias = alias.split(':')[1]
                    if alias in uniprotList:
                        vectorList.append(gene_vec_map[gene[0]])
        patient_vector_list.append([patientPathway[0],vectorList])
        count +=1

    path = os.path.join(config.data_dir, 'KircVectors.txt')
    with open(path, 'w') as f:
        for patient in patient_vector_list:
            l = patient[0]
            for vector in patient[1]:
                l += ";" + ','.join(map(str, vector))
            f.write(l+"\n")
    return patient_vector_list

    print("Done")

#get_n2v_representations()
#findPathwaysForAllPatients()


