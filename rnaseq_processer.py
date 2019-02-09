import numpy as np
import pandas as pd

def process(fn):
	'''
	Inputs
	fn: csv file in which the first row contains gene name, entrez gene id and patient ids.
	    The rest of the rows are the genes and their expressions (binary value) on each patient.

	Outputs
	gene_expressions: a dataframe indicating the over- and under-expressed genes where
		genes with entrez gene id are on rows and patient ids are on columns
		gene_names: a dataframe indicating the name of the genes of entrez gene ids
	'''

	data = pd.read_csv(fn, sep="\t")
	#data = data.set_index(["Gene Name", "Entrez Gene ID"])
	#data = data.set_index("Entrez Gene ID")
	data[["gene_name", "entrez_gene_id"]] = data['#probe'].str.split('|', n=1, expand=True)
	data = data.set_index(["entrez_gene_id"])
	data = data.drop(columns=["#probe"])
	#print(data.head())
	#print(data == "?")

	# IF THERE ARE ADDITIONAL PROBLEMS TO BE CONSIDERED, PREPROCESS THE DATA ACCORDINGLY
	# drop genes with name "?"
	tmp = data.index[data["gene_name"] == "?"]
	data = data.drop(tmp)

	# sort the dataframe based on both gene name and patient id
	data = data.sort_index(axis=1)
	data = data.sort_values("gene_name")

	gene_names = data["gene_name"]
	data = data.drop(columns=["gene_name"])

	#print(data.columns[0:-1])

	# delete patients with non-solid tumor
	tmp = [row[3][0:2] == '01' for row in data.columns.str.split("-").tolist()]
	#print(tmp[0:5])
	ind = [i for i,x in enumerate(tmp) if x==False] # find the indices of non-solid tumors
	#print(data.columns[ind])
	data = data.drop(columns=data.columns[ind]) # drop them from the data frame
	#print(data)

	# drop the genes which are not expressed more than half of the samples
	genes_to_drop = data.index[(data == 0).T.sum().values > (len(data.columns)/2)]
	data = data.drop(genes_to_drop)
	#print(data)
	#print("\n\n")

	# calculate z-scores
	mean_exp = data.mean(axis=1, numeric_only=True)
	#print(mean_exp)
	std_exp = data.std(axis=1, numeric_only=True)
	#print(std_exp)
	z_scores = data.subtract(mean_exp, axis=0)
	#print(z_scores)
	z_scores = z_scores.div(std_exp, axis=0)
	#print(z_scores)
	#print("\n\n")

	# find differentially expressed genes
	# 1 for over-expressed, -1 for under-expressed
	gene_expression = pd.DataFrame(0, index=data.index, columns=data.columns)
	threshold = 1.96 # two standard deviation
	gene_expression[z_scores > threshold] = 1
	#print(gene_expression)
	gene_expression[z_scores < (-1 * threshold)] = -1
	#gene_expression = gene_expression.sort_index(axis=1)
	#print(gene_expression)

	return gene_expression, gene_names
