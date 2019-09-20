import numpy as np
import pandas as pd
import csv
import os
import config

folder_loc = os.path.join(config.data_dir, 'rppa_converter')


def pruneProteins(data):
    tmp = data['#probe'].str.split('-')
    prots = []
    un_proc_prots = []
    for prot in tmp:
        pieces = prot[:-2]

        un_proc_prot_name = '-'.join(pieces)
        tmp = ''.join(pieces)
        protName = ''.join(tmp.split('_'))
        prots.append(protName.upper())
        un_proc_prots.append(un_proc_prot_name)
    return prots, un_proc_prots


def getEntrezData():
    loc = folder_loc
    p2eg_loc = loc+"/gene_uniprot_entrez_kegg.csv"
    add_map_loc = loc+"/external_mapping_of_rppa_proteins.csv"
    ext_eg_id_map_loc = loc+"/external_mapping_of_rppa_proteins_with_entrez_gene_ids.txt"
    manual_map_loc = loc+"/manual_mapping.txt"


    #First mapping
    p2eg = []
    delimiter = ";"
    start_row = 2
    with open(p2eg_loc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                tempList = [row[0]]+[row[3]]
                p2eg.append(tempList)


    #Second mapping
    add_map = []
    delimiter = ","
    start_row = 1

    with open(add_map_loc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                tempList = [row[0]]+[row[3]]
                add_map.append(row)


    #Third mapping
    ext_eg_id_map = []
    delimiter = ","
    start_row = 1

    with open(ext_eg_id_map_loc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                ext_eg_id_map.append(row)



    #Fourth mapping
    manual_map = []
    delimiter = ","
    start_row = 1

    with open(manual_map_loc, "r", encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row-1):
            next(csv_reader)
        for row in csv_reader:
            if len(row) > 1:
                manual_map.append(row)

    return np.array(p2eg), np.array(add_map), np.array(ext_eg_id_map), np.array(manual_map)

def get_or(x):
    if any(t == 1 for t in x):
        return 1
    elif any(t == -1 for t in x):
        return -1
    else:
        return 0

def process_phos(fn):
    '''
        Inputs
            fn: file loc of synapse rppa data

        Outputs
            gene_expressions: a dataframe indicating the over- (1) and under-expressed (-1) genes where
                genes with entrez gene id are on rows and patient ids are on columns
        '''

    data = pd.read_csv(fn, sep="\t")
    # data = data.set_index(["Gene Name", "Entrez Gene ID"])
    # data = data.set_index("Entrez Gene ID")
    drop_ids = []
    for idx, col in enumerate(data.columns):
        if not (idx == 0) and not ('01' in col.split('-')[3]):
            drop_ids.append(idx)

    data.drop(drop_ids)
    proc_prots, unproc_prots = pruneProteins(data)

    # mappings
    p2eg, add_map, ext_eg_id_map, manual_map = getEntrezData()

    eg_ids = np.zeros(len(proc_prots))
    # non_eg_ids = np.zeros(len(procProts))

    proteins = unproc_prots
    patients = data.columns[1:]

    # Try p2eg
    for i in range(len(data)):
        if proc_prots[i] in p2eg[:, 0]:
            idx = np.where(p2eg[:, 0] == proc_prots[i])
            eg_ids[i] = p2eg[:, 1].astype(int)[idx]

    # Try ext_eg_id_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in ext_eg_id_map[:, 0]:
                idx = np.where(ext_eg_id_map[:, 0] == unproc_prots[i])
                eg_ids[i] = ext_eg_id_map[:, 2][idx]

    # Try manual_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in manual_map[:, 0]:
                idx = np.where(manual_map[:, 0] == unproc_prots[i])
                eg_ids[i] = manual_map[:, 2][idx]

    # data.replace(list(data["#probe"]), list(eg_ids))
    data["#probe"] = eg_ids.astype(int)

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

    data = data.sort_index(axis=1)
    data = data.sort_values('#probe', ascending=True)
    #data = data.drop_duplicates(subset='#probe', keep='first')apply(get_or)
    data["#probe"] = np.array(data["#probe"]).astype('str')
    data = data.set_index(["#probe"])

    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    z_scores = z_scores.div(std_exp, axis=0)

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
    threshold = 1.96  # two standard deviation
    gene_expression[z_scores > threshold] = 1
    gene_expression[z_scores < (-1 * threshold)] = -1
    gene_expression = gene_expression.groupby(gene_expression.index).agg(lambda x: get_or(x))
    #gene_expression = gene_expression.sort_index(axis=1)

    return gene_expression

def process_phos_freq(fn):
    '''
        Inputs
            fn: file loc of synapse rppa data

        Outputs
            gene_expressions: a dataframe indicating the over- (1) and under-expressed (-1) genes where
                genes with entrez gene id are on rows and patient ids are on columns
        '''

    data = pd.read_csv(fn, sep="\t")
    # data = data.set_index(["Gene Name", "Entrez Gene ID"])
    # data = data.set_index("Entrez Gene ID")
    drop_ids = []
    for idx, col in enumerate(data.columns):
        if not (idx == 0) and not ('01' in col.split('-')[3]):
            drop_ids.append(idx)

    data.drop(drop_ids)
    proc_prots, unproc_prots = pruneProteins(data)

    # mappings
    p2eg, add_map, ext_eg_id_map, manual_map = getEntrezData()

    eg_ids = np.zeros(len(proc_prots))
    # non_eg_ids = np.zeros(len(procProts))

    proteins = unproc_prots
    patients = data.columns[1:]

    # Try p2eg
    for i in range(len(data)):
        if proc_prots[i] in p2eg[:, 0]:
            idx = np.where(p2eg[:, 0] == proc_prots[i])
            eg_ids[i] = p2eg[:, 1].astype(int)[idx]

    # Try ext_eg_id_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in ext_eg_id_map[:, 0]:
                idx = np.where(ext_eg_id_map[:, 0] == unproc_prots[i])
                eg_ids[i] = ext_eg_id_map[:, 2][idx]

    # Try manual_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in manual_map[:, 0]:
                idx = np.where(manual_map[:, 0] == unproc_prots[i])
                eg_ids[i] = manual_map[:, 2][idx]

    # data.replace(list(data["#probe"]), list(eg_ids))
    data["#probe"] = eg_ids.astype(int)

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

    data = data.sort_index(axis=1)
    data = data.sort_values('#probe', ascending=True)
    #data = data.drop_duplicates(subset='#probe', keep='first')apply(get_or)
    data["#probe"] = np.array(data["#probe"]).astype('str')
    data = data.set_index(["#probe"])
    uniq_probs = np.unique(data.index)
    for idx in uniq_probs:
        values = data.loc[idx].values
        lenn = sum(data.index == idx)
        if lenn > 1:
            divider = values[0]
            new_rows = np.zeros((lenn,len(values[0])))
            new_rows[0] = divider
            for i in range(1,lenn):
                new_rows[i] = np.divide(values[i],divider)
            data.loc[idx] = new_rows
    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    z_scores = z_scores.div(std_exp, axis=0)

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
    threshold = 1.96  # two standard deviation
    gene_expression[z_scores > threshold] = 1
    gene_expression[z_scores < (-1 * threshold)] = -1
    gene_expression = gene_expression.groupby(gene_expression.index).agg(lambda x: get_or(x))
    #gene_expression = gene_expression.sort_index(axis=1)

    return gene_expression


def process(fn):
    '''
    Inputs
        fn: file loc of synapse rppa data

    Outputs
        gene_expressions: a dataframe indicating the over- (1) and under-expressed (-1) genes where
            genes with entrez gene id are on rows and patient ids are on columns
    '''

    data = pd.read_csv(fn, sep="\t")
    # data = data.set_index(["Gene Name", "Entrez Gene ID"])
    # data = data.set_index("Entrez Gene ID")
    drop_ids = []
    for idx, col in enumerate(data.columns):
        if not(idx == 0) and not('01' in col.split('-')[3]):
            drop_ids.append(idx)

    data.drop(drop_ids)
    proc_prots, unproc_prots = pruneProteins(data)

    # mappings
    p2eg, add_map, ext_eg_id_map, manual_map = getEntrezData()

    eg_ids = np.zeros(len(proc_prots))
    # non_eg_ids = np.zeros(len(procProts))


    proteins = unproc_prots
    patients = data.columns[1:]

    # Try p2eg
    for i in range(len(data)):
        if proc_prots[i] in p2eg[:, 0]:
            idx = np.where(p2eg[:, 0] == proc_prots[i])
            eg_ids[i] = p2eg[:, 1].astype(int)[idx]

    # Try ext_eg_id_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in ext_eg_id_map[:,0]:
                idx = np.where(ext_eg_id_map[:, 0] == unproc_prots[i])
                eg_ids[i] = ext_eg_id_map[:, 2][idx]

    # Try manual_map
    for i in range(len(unproc_prots)):
        if eg_ids[i] == 0:
            if unproc_prots[i] in manual_map[:, 0]:
                idx = np.where(manual_map[:, 0] == unproc_prots[i])
                eg_ids[i] = manual_map[:, 2][idx]

    # data.replace(list(data["#probe"]), list(eg_ids))
    data["#probe"] = eg_ids.astype(int)

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns)/2)]
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

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
    threshold = 1.96  # two standard deviation
    gene_expression[z_scores > threshold] = 1
    gene_expression[z_scores < (-1 * threshold)] = -1

    return gene_expression

#process_phos("../data/kirc_data/kirc_rppa_data")