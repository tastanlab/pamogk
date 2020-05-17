import csv

import numpy as np
import pandas as pd

from .. import config

RPPA_DATA_DIR = config.DATA_DIR / 'rppa'


def prune_proteins(data):
    tmp = data['#probe'].str.split('-')
    prots = []
    un_proc_prots = []
    for prot in tmp:
        pieces = prot[:-2]

        un_proc_prot_name = '-'.join(pieces)
        tmp = ''.join(pieces)
        prot_name = ''.join(tmp.split('_'))
        prots.append(prot_name.upper())
        un_proc_prots.append(un_proc_prot_name)
    return prots, un_proc_prots


def read_csv(path, has_header=False, delimiter=',', mapper=lambda r: r):
    path = RPPA_DATA_DIR / path
    with open(path, 'r', encoding='utf8', errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        if has_header:
            next(csv_reader)
        return [mapper(r) for r in csv_reader if len(r) > 1]


def get_entrez_data():
    # First mapping
    p2eg = read_csv('gene_uniprot_entrez_kegg.csv', has_header=True, delimiter=';', mapper=lambda r: [r[0], r[3]])

    # Second mapping
    add_map = read_csv('external_mapping_of_rppa_proteins.csv')

    # Third mapping
    ext_eg_id_map = read_csv('external_mapping_of_rppa_proteins_with_entrez_gene_ids.txt')

    # Fourth mapping
    manual_map = read_csv('manual_mapping.txt')

    return np.array(p2eg), np.array(add_map), np.array(ext_eg_id_map), np.array(manual_map)


def process(filename):
    """
    :param filename: file loc of synapse rppa data
    :returns: a dataframe indicating the over- (1) and under-expressed (-1) genes where
            genes with entrez gene id are on rows and patient ids are on columns
    """

    data = pd.read_csv(config.get_safe_data_file(filename), sep='\t')
    # data = data.set_index(['Gene Name', 'Entrez Gene ID'])
    # data = data.set_index('Entrez Gene ID')
    drop_cols = []
    for col in data.columns:
        if not (col == '#probe') and not ('01' in col.split('-')[3]):
            drop_cols.append(col)

    data.drop(columns=drop_cols)
    proc_prots, unproc_prots = prune_proteins(data)

    # mappings
    p2eg, add_map, ext_eg_id_map, manual_map = get_entrez_data()

    eg_ids = np.zeros(len(proc_prots))
    # non_eg_ids = np.zeros(len(procProts))

    # Try p2eg
    for i in range(len(data)):
        if proc_prots[i] in p2eg[:, 0]:
            idx = np.where(p2eg[:, 0] == proc_prots[i])
            eg_ids[i] = p2eg[:, 1].astype(int)[idx]

    def map_eg_id(src):
        for i in range(len(unproc_prots)):
            if eg_ids[i] == 0:
                if unproc_prots[i] in src[:, 0]:
                    idx = np.where(src[:, 0] == unproc_prots[i])
                    eg_ids[i] = src[:, 2][idx]

    # Try ext_eg_id_map
    map_eg_id(ext_eg_id_map)

    # Try manual_map
    map_eg_id(manual_map)

    # data.replace(list(data['#probe']), list(eg_ids))
    data['#probe'] = eg_ids.astype(int)

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

    data = data.sort_index(axis=1)
    data = data.sort_values('#probe', ascending=True)
    data = data.drop_duplicates(subset='#probe', keep='first')
    data['#probe'] = np.array(data['#probe']).astype('str')
    data = data.set_index(['#probe'])

    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    z_scores = z_scores.div(std_exp, axis=0)

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
    threshold = 1.96  # t   wo standard deviation
    gene_expression[z_scores > threshold] = 1
    gene_expression[z_scores < (-1 * threshold)] = -1

    return gene_expression


def process_cont(filename):
    """
    :param filename: file loc of synapse rppa data
    :return: a dataframe indicating the normalized values of expressions of genes where
                genes with entrez gene id are on rows and patient ids are on columns
    """
    data = pd.read_csv(config.get_safe_data_file(filename), sep='\t')
    # data = data.set_index(['Gene Name', 'Entrez Gene ID'])
    # data = data.set_index('Entrez Gene ID')
    drop_cols = []
    for col in data.columns:
        if not (col == '#probe') and not ('01' in col.split('-')[3]):
            drop_cols.append(col)

    data.drop(columns=drop_cols)
    proc_prots, unproc_prots = prune_proteins(data)

    # mappings
    p2eg, add_map, ext_eg_id_map, manual_map = get_entrez_data()

    eg_ids = np.zeros(len(proc_prots))
    # non_eg_ids = np.zeros(len(procProts))

    proteins = unproc_prots
    patients = data.columns[1:]

    # Try p2eg
    for i in range(len(data)):
        if proc_prots[i] in p2eg[:, 0]:
            idx = np.where(p2eg[:, 0] == proc_prots[i])
            eg_ids[i] = p2eg[:, 1].astype(int)[idx]

    def map_eg_id(src):
        for i in range(len(unproc_prots)):
            if eg_ids[i] == 0:
                if unproc_prots[i] in src[:, 0]:
                    idx = np.where(src[:, 0] == unproc_prots[i])
                    eg_ids[i] = src[:, 2][idx]

    # Try ext_eg_id_map
    map_eg_id(ext_eg_id_map)

    # Try manual_map
    map_eg_id(manual_map)

    # data.replace(list(data["#probe"]), list(eg_ids))
    data["#probe"] = eg_ids.astype(int)

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

    data = data.sort_index(axis=1)
    data = data.sort_values('#probe', ascending=True)
    data = data.drop_duplicates(subset='#probe', keep='first')
    data["#probe"] = np.array(data["#probe"]).astype('str')
    data = data.set_index(["#probe"])

    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    z_scores = z_scores.div(std_exp, axis=0)

    # find normalized expressions of genes
    gene_expression = pd.DataFrame(z_scores, index=data.index, columns=data.columns)

    return gene_expression
