import pandas as pd

from .. import config


def process(filename):
    """
    :param filename: tab separated file in which the first row contains gene name/entrez gene id combined
            and patient ids. The rest of the rows are the genes and their expressions on each patient.

    :returns
        gene_expressions: a dataframe indicating the over- (1) and under-expressed (-1) genes where
            genes with entrez gene id are on rows and patient ids are on columns
        gene_name_map: a dataframe indicating the name of the genes of entrez gene ids
    """

    fpath = config.get_safe_data_file(filename)
    data = pd.read_csv(fpath, sep='\t')
    # data = data.set_index(['Gene Name', 'Entrez Gene ID'])
    # data = data.set_index('Entrez Gene ID')
    data[['gene_name', 'entrez_gene_id']] = data['#probe'].str.split('|', n=1, expand=True)
    data = data.set_index(['entrez_gene_id'])
    data = data.drop(columns=['#probe'])

    # IF THERE ARE ADDITIONAL PROBLEMS TO BE CONSIDERED, PREPROCESS THE DATA ACCORDINGLY
    # drop genes with name '?'
    tmp = data.index[data['gene_name'] == '?']
    data = data.drop(tmp)

    # sort the dataframe based on both gene name and patient id
    data = data.sort_index(axis=1)
    data = data.sort_values('gene_name')

    gene_name_map = data['gene_name']
    data = data.drop(columns=['gene_name'])

    # delete patients with non-solid tumor
    tmp = [row[3][0:2] == '01' for row in data.columns.str.split('-').tolist()]
    ind = [i for i, x in enumerate(tmp) if x == False]  # find the indices of non-solid tumors
    data = data.drop(columns=data.columns[ind])  # drop them from the data frame

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

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

    return gene_expression, gene_name_map


def process_cont(fn):
    """
    Inputs
        fn: tab separated file in which the first row contains gene name/entrez gene id combined
        and patient ids. The rest of the rows are the genes and their expressions on each patient.

    Outputs
        gene_expressions: a dataframe indicating the normalized values of expressions of genes where
                genes with entrez gene id are on rows and patient ids are on columns
        gene_name_map: a dataframe indicating the name of the genes of entrez gene ids
    """

    data = pd.read_csv(fn, sep='\t')
    # data = data.set_index(['Gene Name', 'Entrez Gene ID'])
    # data = data.set_index('Entrez Gene ID')
    data[['gene_name', 'entrez_gene_id']] = data['#probe'].str.split('|', n=1, expand=True)
    data = data.set_index(['entrez_gene_id'])
    data = data.drop(columns=['#probe'])

    # IF THERE ARE ADDITIONAL PROBLEMS TO BE CONSIDERED, PREPROCESS THE DATA ACCORDINGLY
    # drop genes with name '?'
    tmp = data.index[data['gene_name'] == '?']
    data = data.drop(tmp)

    # sort the dataframe based on both gene name and patient id
    data = data.sort_index(axis=1)
    data = data.sort_values('gene_name')

    gene_name_map = data['gene_name']
    data = data.drop(columns=['gene_name'])

    # delete patients with non-solid tumor
    tmp = [row[3][0:2] == '01' for row in data.columns.str.split('-').tolist()]
    ind = [i for i, x in enumerate(tmp) if x == False]  # find the indices of non-solid tumors
    data = data.drop(columns=data.columns[ind])  # drop them from the data frame

    # drop the genes which are not expressed more than half of the samples
    genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns) / 2)]
    data = data.drop(genes_to_drop)

    # calculate z-scores
    mean_exp = data.mean(axis=1, numeric_only=True)
    std_exp = data.std(axis=1, numeric_only=True)
    z_scores = data.subtract(mean_exp, axis=0)
    # print(z_scores)
    z_scores = z_scores.div(std_exp, axis=0)

    # find differentially expressed genes
    # 1 for over-expressed, -1 for under-expressed
    gene_expression = pd.DataFrame(z_scores, index=data.index, columns=data.columns)
    # threshold = 1.96 # two standard deviation
    # gene_expression[z_scores > threshold] = 1
    # gene_expression[z_scores < (-1 * threshold)] = -1

    return gene_expression, gene_name_map
