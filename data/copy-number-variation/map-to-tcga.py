#!/usr/bin/env python
import csv
import requests

pat_ids = []
with open('KIRC.focal_score_by_genes.tsv') as f:
    pat_ids = f.readline().strip().split('\t')[3:]

num_pats = len(pat_ids)
print('Found patient ids of size:', num_pats)
with open('tcga-mapping.csv-2', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['ensembl_gene_id', 'entrez_gene_id'])
    i = 1
    for pat_id in pat_ids:
        res = requests.get(f"https://api.gdc.cancer.gov/v0/all?query={pat_id}&size=5").json()
        hits = res['data']['query']['hits']
        ent_id = hits[0]['submitter_id'] if len(hits) > 0 else None
        writer.writerow([pat_id, ent_id])
        print(f'\tFound match {i:3}/{num_pats} {pat_id} --> {ent_id}\r', end='')
        i += 1
    print('\n')
