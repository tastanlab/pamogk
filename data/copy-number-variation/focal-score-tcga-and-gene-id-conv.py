#!/usr/bin/env python3
import csv
import pdb

import sys

gene_type = 'uniprot'
if len(sys.argv) > 1 and sys.argv[1] == 'entrez':
    gene_type = 'entrez'

print('Loading patient mapping from case UUID to tcga')
pat_mapping = {}
with open('tcga-mapping.csv') as f:
    for row in f.readlines():
        uuid, pid = row.strip().split(',')
        if uuid in pat_mapping:
            print(f'Duplicate mapping {uuid} --> {pat_mapping[uuid]}, {pid}')
            raise
        pat_mapping[uuid] = pid
print('Loaded unique patient mappings:', len(pat_mapping))
print()

print(f'Loading gene mapping from ensembl to {gene_type}')
gene_mapping = {}
with open(f'ens-to-{gene_type}-mapping.csv') as f:
    for row in f.readlines():
        row_spt = row.strip().split(',')
        if len(row_spt) < 2: continue
        uuid, gid = row_spt
        if uuid not in gene_mapping:
            gene_mapping[uuid] = [gid]
        else:
            gene_mapping[uuid].append(gid)
    print()
print('Loaded unique gene mappings:', len(gene_mapping))
print()

pat_headers = None
skipped = 0
added = 0
with open('KIRC.focal_score_by_genes.tsv') as f, open(f'KIRC.focal_score_by_genes-tcga_{gene_type}.csv', 'w') as q:
    csvReader = csv.DictReader(f, delimiter='\t')
    csvWriter = csv.writer(q)
    for row in csvReader:
        if pat_headers is None:
            pat_headers = [k for k in row.keys() if k in pat_mapping]
            headers = ['Entrez Gene ID'] + [pat_mapping[k] for k in pat_headers]
            csvWriter.writerow(headers)
        gid = row['Gene Symbol'].split('.')[0]
        if gid not in gene_mapping:
            print(f'Skipping gene with id {gid} not found in mapping\r', end='')
            skipped += 1
            continue
        for eid in gene_mapping[gid]:
            added += 1
            csvWriter.writerow([eid] + [row[k] for k in pat_headers])
    print()

print(f'Added: {added}\nSkipped: {skipped}')
