#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
from pathlib import Path

from .. import config
from ..lib.sutils import ensure_file_dir

KIRC_DATA_DIR = config.DATA_DIR / 'kirc_data'

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--somatic-data', '-r', metavar='file-path', dest='somatic_data', type=Path, help='Somatic Data',
                    default=KIRC_DATA_DIR / 'kirc_somatic_mutation_data.csv')


def read_csv(filename, delimiter=','):
    with open(KIRC_DATA_DIR / filename, encoding='utf8') as f:
        return [r for r in csv.DictReader(f, delimiter=delimiter)]


def transpose_data(data, id_col, ignored_cols=None):
    """
    transposes given data from patient columns to patient rows
    data: list
    id_col: str
    ignored_cols: list
    :return: list
    """
    if ignored_cols is None:
        ignored_cols = []
    ignored_cols.append(id_col)  # add id col to ignores
    tm = {}
    pat_ids = set(data[0].keys()).difference(ignored_cols)
    for p in pat_ids:  # initialize patients
        tm[p] = {'Patient ID': p}
    for r in data:
        col_id = r[id_col]
        for p in pat_ids:
            tm[p][col_id] = r[p]
    return list(tm.values())


def add_data_to_patient(key, data, pat_data=None):
    """
    adds kernel data to data map

    :param key:
    :param data:
    :param pat_data:
    :return:
    """
    if pat_data is None:
        pat_data = {}
    for r in data:
        # patient ids for rppa have some extra fields
        pat_id = r['Patient ID'] = '-'.join(r['Patient ID'].split('-')[:3])
        if pat_id not in pat_data:
            pat_data[pat_id] = {}
        pat_data[pat_id][key] = r
    return pat_data


def write_csv(data, key, pat_ids, use_index=False):
    filepath = KIRC_DATA_DIR / 'rmkl' / f'{key}.csv'
    print('Writing patient data type', key, 'to path', filepath)
    ensure_file_dir(filepath)
    with open(filepath, 'w') as f:
        key_data = [r[key] for r in data]
        if use_index:
            for r in key_data:
                r['Patient ID'] = pat_ids.index(r['Patient ID'])
        # make sure patient id is the first column
        fieldnames = ['Patient ID'] + list(k for k in key_data[0].keys() if k not in ['Patient ID'])
        csv_writer = csv.DictWriter(f, fieldnames)
        csv_writer.writeheader()
        csv_writer.writerows(key_data)


def write_pat_id_map(full_pat_ids):
    filepath = KIRC_DATA_DIR / 'rmkl' / 'pat-id-map.csv'
    print('Writing patient id map to path', filepath)
    ensure_file_dir(filepath)
    with open(filepath, 'w') as f:
        csv_writer = csv.DictWriter(f, ['Patient ID', 'index'])
        csv_writer.writeheader()
        for i, pat_id in enumerate(full_pat_ids):
            csv_writer.writerow({'Patient ID': pat_id, 'index': i})


def main():
    clinical_data = read_csv('kirc_clinical_data.csv')
    # somatic_mut_data = read_csv('kirc_somatic_mutation_data.csv')
    rnaseq_data = transpose_data(read_csv('kirc_rna_seq_expression_data.csv'), 'Entrez Gene ID', ['Gene Name'])
    rppa_data = transpose_data(read_csv('kirc_rppa_data', delimiter='\t'), '#probe')

    data = add_data_to_patient('clinical', clinical_data)
    # data = add_data_to_patient('somatic', somatic_mut_data, data)
    data = add_data_to_patient('rnaseq', rnaseq_data, data)
    data = add_data_to_patient('rppa', rppa_data, data)

    print('Total patients:', len(data))

    # get patients that have all 4 data
    full_pats_ids = [k for k, v in data.items() if len(v) == 3]
    full_pats = [data[pat_id] for pat_id in full_pats_ids]
    print('Patients with full data:', len(full_pats))

    write_pat_id_map(full_pats_ids)
    write_csv(full_pats, 'clinical', full_pats_ids)
    # write_csv(full_pats, 'somatic', full_pats_ids)
    write_csv(full_pats, 'rnaseq', full_pats_ids)
    write_csv(full_pats, 'rppa', full_pats_ids)

    # illimuna_data = read_csv('unc.edu_KIRC_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor.txt')


if __name__ == '__main__':
    main()
