import csv
import pdb
import shutil
from pathlib import Path
import sys

import biomart

server = biomart.BiomartServer("http://useast.ensembl.org/biomart")
server.verbose = False
print('BioMart server alive:', server.is_alive)
# accessing first db than dataset makes it faster otherwise it fetches all datasets per db
db = server.databases['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
print('Loaded DB')
ens_genes = []  # ['ENSG00000008128']

args = sys.argv[1:]
target_gene_id = 'entrezgene_id'
MAPPING_FILE = Path('ens-to-entrez-mapping.csv')

if len(args) > 0 and args[0] == 'uniprot':
    target_gene_id = 'uniprot_gn_id'
    MAPPING_FILE = Path('ens-to-uniprot-mapping.csv')

print('Matching with target gene id of type:', target_gene_id)

with open('KIRC.focal_score_by_genes.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        row['Gene'] = row['Gene Symbol'].split('.')[0]
        ens_genes.append(row)
print('Looking up for gene set of size:', len(ens_genes))

pdb.set_trace()

def get_query(ens_genes):
    return {
        'filters': {
            'ensembl_gene_id': [e['Gene'] for e in ens_genes],
        },
        'attributes': ['ensembl_gene_id', target_gene_id],
    }


if MAPPING_FILE.exists():
    now = datetime.now().strftime('%Y%m%d%H%M%S')
    shutil.copyfile(MAPPING_FILE, MAPPING_FILE.with_name(MAPPING_FILE.name + '-' + now))

with open(MAPPING_FILE, 'w') as f:
    writer = csv.writer(f)
    n = 100
    for i in range(0, len(ens_genes), n):
        print(f'\tFetching in range ({i}, {i + n})')
        is_first = 1 if i == 0 else 0
        query = get_query(ens_genes[i:i + n])
        response = db.search(query, header=is_first)
        writer.writerows(line.decode('utf-8').strip().split('\t') for line in response.iter_lines())
        f.flush()
