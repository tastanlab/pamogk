"""
    This file is to process clincal file
    It is enough to call processOneClinicalSomatic with file location as argument to function:
        import synapse_som_processor as ssp
        ssp.processOneCancerSomatic('folder/fileName.maf')
        This returns 2d numpy array with columns PatientId, AliveOrDead(0 or 1),
         days_to_last_followup or days_to_death(alive or dead), Stage

    To write the report to a file:
        import processors.synapse_clinical_processor as clinic
        output = clinical.processOneClinicalSomatic('folder/fileName.maf')
        ssp.writeToFile(output, fileLocation)



"""

import csv

import numpy as np

from .. import config

CANCER_TYPES = ['BLCA', 'BRCA', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'READ', 'UCEC']


def process_one_clinical_somatic(filepath, delimiter='\t', start_row=1):
    data_array = []

    with open(filepath, 'r', encoding='utf8', errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(start_row):
            header = next(csv_reader)
        for row in csv_reader:
            data_array.append(row)

    if 'tumor_stage' in header:
        ts_key = 'tumor_stage'
    elif 'ajcc_neoplasm_disease_stage' in header:
        ts_key = 'ajcc_neoplasm_disease_stage'
    elif 'gynecologic_tumor_grouping_figo_stage' in header:
        ts_key = 'gynecologic_tumor_grouping_figo_stage'
    else:
        ts_key = '#'
    data_indexes = [header.index('#'), header.index('vital_status'), header.index('days_to_death'),
                    header.index('days_to_last_followup'), header.index(ts_key)]

    '''
    for idx, col in enumerate(header):
        if col == '#':
            data_indexes[0] = idx
        if col == 'tumor_stage':
            data_indexes[4] = idx
        if col == 'vital_status':
            data_indexes[1] = idx
        if col == 'days_to_death':
            data_indexes[2] = idx
        if col == 'days_to_last_followup':
            data_indexes[3] = idx

    '''

    data_array_processed = []
    for row in data_array:
        tmp = row[data_indexes[0]]
        splitted = tmp.split('-')
        if (len(splitted) > 3) and ('01' in splitted[3]):
            new_val = '-'
            new_val = new_val.join(splitted[0:3])
            if row[data_indexes[1]] == 'LIVING':
                temp_days = row[data_indexes[3]]
                row3_value = 0
            else:
                temp_days = row[data_indexes[2]]
                row3_value = 1
            stage = 'NA'
            if row[data_indexes[4]] != 'NA':
                checker = row[data_indexes[4]]
                if ('iv' in checker) or ('IV' in checker):
                    stage = 4
                elif ('iii' in checker) or ('III' in checker):
                    stage = 3
                elif ('ii' in checker) or ('II' in checker):
                    stage = 2
                elif ('i' in checker) or ('I' in checker):
                    stage = 1
                '''
                listFrom = ['Stage I', 'Stage II', 'Stage III', 'Stage IV']
                if row[data_indexes[4]] in listFrom:
                    indforTo = listFrom.index(row[data_indexes[4]])
                    listTo = [1, 2, 3, 4]
                    stage = listTo[indforTo]
                '''
            temp_list = [new_val] + [row3_value] + [temp_days] + [stage]
            data_array_processed.append(temp_list)

    return np.array(data_array_processed)


# 'kirc/data/synapse_kirc_som_data.csv'
def write_to_file(report, filepath):
    with open(filepath, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['Patient ID', 'Status', 'Days', 'Stage'])
        csv_writer.writerows(report)


def print_report(report):
    if len(report) != 0:
        print(len(report))
        print(len(list(set(report[:, 0]))))
        alive_count = np.count_nonzero(report[:, 1] == '0')
        dead_count = np.count_nonzero(report[:, 1] == '1')
        stage_counts = [np.count_nonzero(report[:, 3] == str(i + 1)) for i in range(4)]
        print(alive_count)
        print(dead_count)
        for s in stage_counts: print(s)


def report_all_cancer_types(dir_path):
    for ct in CANCER_TYPES:
        filepath = dir_path / ct / 'clinical'
        report = process_one_clinical_somatic(filepath)
        print(ct + '# of row, # of Unique Gene, #of Unique Patient')
        print_report(report)


CANCER_TYPE = 'OV'

FILEPATH = config.ROOT_DIR / CANCER_TYPE / 'clinical'
rep = process_one_clinical_somatic(FILEPATH)
write_to_file(rep, config.DATA_DIR / 'ov_data' / 'ov_clinical_data.csv')

# report_all_cancer_types(DATA_DIR)
