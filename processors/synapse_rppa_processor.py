import numpy as np
import pandas as pd
import csv

folderLoc = "/home/yitepeli/ForExp/rppa/"
def pruneProteins(data):
    tmp = data['#probe'].str.split('-')
    prots = []
    unProcProts = []
    for prot in tmp:
        pieces = prot[:-2]

        unProcProtName = '-'.join(pieces)
        tmp = ''.join(pieces)
        protName = ''.join(tmp.split('_'))
        prots.append(protName.upper())
        unProcProts.append(unProcProtName)
    return prots,unProcProts

def getEntrezData():
    loc = folderLoc
    p2egLoc = loc+"gene_uniprot_entrez_kegg.csv"
    add_mapLoc = loc+"external_mapping_of_rppa_proteins.csv"
    ext_eg_id_mapLoc = loc+"external_mapping_of_rppa_proteins_with_entrez_gene_ids.txt"
    manual_mapLoc = loc+"manual_mapping.txt"


    #First mapping
    p2eg = []
    delimiter = ";"
    startRow = 2
    with open(p2egLoc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                tempList = [row[0]]+[row[3]]
                p2eg.append(tempList)


    #Second mapping
    add_map = []
    delimiter = ","
    startRow = 1

    with open(add_mapLoc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                tempList = [row[0]]+[row[3]]
                add_map.append(row)


    #Third mapping
    ext_eg_id_map = []
    delimiter = ","
    startRow = 1

    with open(ext_eg_id_mapLoc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                ext_eg_id_map.append(row)



    #Fourth mapping
    manual_map = []
    delimiter = ","
    startRow = 1

    with open(manual_mapLoc, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow-1):
            next(csv_reader)
        for row in csv_reader:
            if(len(row)>1):
                manual_map.append(row)

    return np.array(p2eg), np.array(add_map), np.array(ext_eg_id_map), np.array(manual_map)

def process(fn):
    '''
    Inputs
        fn: file loc of synapse rppa data

    Outputs
        gene_expressions: a dataframe indicating the over- (1) and under-expressed (-1) genes where
            genes with entrez gene id are on rows and patient ids are on columns
    '''

    data = pd.read_csv(fn, sep="\t")
    #data = data.set_index(["Gene Name", "Entrez Gene ID"])
    #data = data.set_index("Entrez Gene ID")
    dropIds = []
    for idx,col in enumerate(data.columns):
        if not(idx==0) and not('01' in col.split('-')[3]) :
            dropIds.append(idx)

    data.drop(dropIds)
    procProts, unProcProts = pruneProteins(data)

    #mappings
    p2eg, add_map, ext_eg_id_map, manual_map = getEntrezData()

    eg_ids = np.zeros(len(procProts))
    #non_eg_ids = np.zeros(len(procProts))


    proteins = unProcProts
    patients = data.columns[1:]

    #Try p2eg
    for i in range(len(data)):
        if procProts[i] in p2eg[:,0]:
            idx = np.where(p2eg[:,0] == procProts[i])
            eg_ids[i] = p2eg[:,1].astype(int)[idx]

    #Try ext_eg_id_map
    for i in range(len(unProcProts)):
        if eg_ids[i]==0:
            if unProcProts[i] in ext_eg_id_map[:,0]:
                idx = np.where(ext_eg_id_map[:, 0] == unProcProts[i])
                eg_ids[i] = ext_eg_id_map[:, 2][idx]

    #Try manual_map
    for i in range(len(unProcProts)):
        if eg_ids[i]==0:
            if unProcProts[i] in manual_map[:,0]:
                idx = np.where(manual_map[:, 0] == unProcProts[i])
                eg_ids[i] = manual_map[:, 2][idx]

    #data.replace(list(data["#probe"]), list(eg_ids))
    data["#probe"] = eg_ids.astype(int)


    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns)/2)]
    data = data.drop(genes_to_drop)


    data = data.sort_values('#probe', ascending=True)
    data = data.drop_duplicates(subset='#probe', keep='first')
    data = data.set_index(["#probe"])

    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    z_scores = z_scores.div(std_exp, axis=0)

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
    threshold = 1.96 # two standard deviation
    gene_expression[z_scores > threshold] = 1
    gene_expression[z_scores < (-1 * threshold)] = -1

    return gene_expression
#exp = process('/home/yitepeli/ForExp/KIRC/rppa')
#print()