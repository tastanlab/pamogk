# -*- coding: utf-8 -*-
import os
from csv import DictReader
import clinical_config

def filter_ov(r):
    return (r['vital_status'] in ['DECEASED', 'LIVING'] and # vital_status is either living or deceased
        len(r['patient_id']) > 0 and # have patient_id
        r['sample_type'] == 'Primary Tumor' and # primary tumor info
        (len(r['days_to_death']) > 0 or len(r['days_to_last_followup']) > 0) and # one of days_to_death & days_to_last_followup exists
        (len(r['days_to_death']) > 0 or r['vital_status'] != 'DECEASED') and # if deceased then days_to_death should exists
        (len(r['days_to_last_followup']) > 0 or r['vital_status'] != 'LIVING') and # if living then days_to_last_followup should exists
        len(r['tumor_stage'])) # stage

def read_ov():
    fpath = os.path.join(clinical_config.data_root, 'ovarian_files/tcga_OV_clinical.aliquot.whitelist_tumor.txt')
    if not os.path.exists(fpath):
        raise Error('Could not find clinical data on path:' + fpath)
    print('Reading ov data from:', fpath)
    raw_data = []
    with open(fpath) as f:
        csvReader = DictReader(f, delimiter='\t')
        # filtering
        raw_data = [r for r in csvReader]
        tmp_data = [r for r in raw_data if filter_ov(r)]
        print('After filtering rows:', len(raw_data), len(tmp_data))
        raw_data = tmp_data
    # post processing
    stages = {'IA': 1, 'IB': 1, 'IC': 1, 'IIA': 2, 'IIB': 2, 'IIC': 2, 'IIIA': 3, 'IIIB': 3, 'IIIC': 3, 'IV': 4};
    for r in raw_data: r['tumor_stage'] = stages[r['tumor_stage']]
    # extract unique patients and their data
    patients = set()
    tmp_data = []
    for r in raw_data:
        pid = r['patient_id']
        if not pid in patients:
            tmp_data.append(r)
            patients.add(pid)
    print('After unique rows:', len(raw_data), len(tmp_data))
    raw_data = tmp_data

    # create data from patient id - overall survival - vital status
    data = []
    for r in raw_data:
        nr = { 'patient_id': r['patient_id'] }
        if r['vital_status'] == 'DECEASED':
            nr['deceased'] = True
            nr['overall_survival'] = r['days_to_death']
        else:
            nr['deceased'] = False
            nr['overall_survival'] = r['days_to_last_followup']
        data.append(nr)
    return data, raw_data


