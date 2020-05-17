import argparse
import csv
from pathlib import Path

import numpy as np

from .. import config
from ..lib.sutils import ensure_file_dir

'''
    This file is to process som file which has type of .maf
    It is enough to call processOneCancerSomatic with file location as argument to function:
        import synapse_som_processor as ssp
        ssp.processOneCancerSomatic('folder/fileName.maf')
        This returns 2d numpy array with columns GeneName, EntrezId, Patient Id

    To write the report to a file:
        import synapse_som_processor as ssp
        output = ssp.processOneCancerSomatic('folder/fileName.maf')
        ssp.writeToFile(output, fileLocation)
'''

# CancerTypes = ['BLCA','BRCA','COAD','GBM','HNSC','KIRC','LAML','LUAD','LUSC','OV','READ','UCEC']
CANCER_TYPES = ['BLCA', 'COAD', 'GBM', 'HNSC', 'LAML', 'LUAD', 'LUSC', 'OV', 'READ', 'UCEC']
parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--somatic-data', '-r', metavar='file-path', dest='somatic_data', type=Path, help='Somatic Data',
                    default=config.DATA_DIR / 'kirc_data/kirc_somatic_mutation_data.csv')

args = parser.parse_args()


def process_one_cancer_somatic(filepath, delimiter='\t', start_row=1):
    data_array = []

    with open(filepath, 'r', encoding='utf8', errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row):
            next(csv_reader)
        for row in csv_reader:
            if len(row) > 1:
                temp_list = [row[0]] + [row[1]] + [row[15]]
                data_array.append(temp_list)

    output = np.array(data_array)
    prune_ind = np.zeros(len(output))
    prune_list = []
    for idx, data in enumerate(output):
        tmp = data[2]
        splitted = tmp.split('-')
        if '01' in splitted[3]:
            prune_ind[idx] = 1
            data[2] = '-'.join(splitted[0:3])
        else:
            prune_list.append(idx)

    output = np.delete(output, prune_list, axis=0)

    return np.sort(output, axis=0)


# 'kirc/data/synapse_kirc_som_data.csv'
def write_to_file(report, filepath):
    with open(filepath, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['Gene Name', 'Entrez Gene ID', 'Patient ID'])
        csv_writer.writerows(report)


def print_report(report):
    print(len(report))
    print(len(list(set(report[:, 0]))))
    print(len(list(set(report[:, 2]))))


def report_all_cancer_types(data_dir):
    for cancer in CANCER_TYPES:
        filepath = data_dir / cancer / 'som.maf'
        report = process_one_cancer_somatic(filepath)
        print(cancer + '# of row, # of Unique Gene, #of Unique Patient')
        print_report(report)


def read_processed_data():
    # Real Data #
    # process RNA-seq expression data
    patients = {}
    with open(args.somatic_data) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pat_id = row['Patient ID']
            ent_id = row['Entrez Gene ID']
            if pat_id not in patients:
                patients[pat_id] = {ent_id}
            else:
                patients[pat_id].add(ent_id)

    return patients

def process_and_save_cancer(cancer_type='BRCA', data_dir=Path('/home/yitepeli/ForExp/')):
    ct_lower = cancer_type.lower()
    outfile_path = data_dir / f'{ct_lower}_data' / f'{ct_lower}_somatic_mutation_data.csv'
    ensure_file_dir(outfile_path)
    filepath = data_dir / cancer_type / 'som.maf'
    rep = process_one_cancer_somatic(filepath)
    write_to_file(rep, outfile_path)


def process_and_save_all():
    for ct in CANCER_TYPES:
        process_and_save_cancer(ct)


# patient_dict = read_processed_data()
# draw_hist(patient_dict)
# report_all_cancer_types()
process_and_save_all()
print()
