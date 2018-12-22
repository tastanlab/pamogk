#!/usr/bin/env python3
import requests
import config
import numpy as np
import os
import pdb
from pathway_reader import cx_pathway_reader
import json

def get_all_nodes(all_pws_map):
    return np.hstack([[g.nodes[n] for n in g.nodes()] for g in all_pws_map.values()])

def get_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.hstack([n['alias'] for n in all_nodes if 'alias' in n])

def get_non_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.hstack([n['n'] for n in all_nodes if 'alias' not in n])

def fetch_uniprot_mapping(alias_list, fr='ACC,ID', to='P_ENTREZGENEID'):
    stringList = ",".join(alias_list)
    print(stringList)
    r = requests.post('https://www.uniprot.org/uploadlists/',
        headers = {
            'authority': 'www.uniprot.org',
            'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36',
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
            'upgrade-insecure-requests': '1',
            'referer': 'https://www.uniprot.org/uploadlists/',
        },
        data = {
            'landingPage': 'false',
            'jobId': '',
            # 'uploadQuery': ' '.join(alias_list),
            'uploadQuery': stringList,
            'url2': '',
            'from': fr,
            'to': to,
            'taxon': '',
        })
    MAPPING_ID = None
    if len(r.history) > 0: MAPPING_ID = r.history[0].headers['location'].split('/')[-1]
    if MAPPING_ID is None:
        raise Exception('Could not find mapping for aliases')
    r1 = requests.get('https://www.uniprot.org/mapping/%s.tab' % MAPPING_ID)
    print(r1.text)
    r2 = requests.get('https://www.uniprot.org/mapping/%s.not' % MAPPING_ID)
    print(r2.text)
    return r1.text, r2.text
    # pdb.set_trace()

def fetch_uniprot_to_entrez(alias_list):
    uniprot_alias_list = [a[1] for a in [a.split('uniprot knowledgebase:')for a in alias_list] if len(a) > 1]
    return fetch_uniprot_mapping(uniprot_alias_list)

def get_uniprot_to_entrez_map(all_pws_map = None):
    if all_pws_map is None:
        all_pws_map = cx_pathway_reader.read_pathways()
    all_alias = get_alias_set(all_pws_map)
    data = fetch_uniprot_to_entrez(all_alias)
    mappedOnes = data[0]
    notMappedOnes = data[1]
    jsonData = {}
    jsonData['mapped'] = []
    jsonData['notMapped'] = []


    mappedList = mappedOnes.split("\n")
    for i in range(1,len(mappedList)-1):
        jsonData['mapped'].append({
            'from': mappedList[i].split('\t')[0],
            'to':mappedList[i].split('\t')[1]
        })

    notMappedList = notMappedOnes.split("\n")
    for i in range(1,len(notMappedList)-1):
        jsonData['notMapped'].append({
            'from': notMappedList[i],
            'to':''
        })

    path = os.path.join(config.data_dir, 'gene-map-uniprot-entrez.json')
    with open(path, 'w') as f:
        json.dump(jsonData, f)

def jsonToDict():
    path = os.path.join(config.data_dir, 'gene-map-uniprot-entrez.json')
    with open(path,'r') as f:
        data = json.load(f)

    mappedData = data['mapped']
    notMappedData = data['notMapped']

    uniProtToEntrez = {}
    entrezToUniProt = {}

    for line in mappedData:
        if line['from'] in uniProtToEntrez:
            uniProtToEntrez[line['from']].append(line['to'])
        else:
            uniProtToEntrez[line['from']] = [line['to']]

        if line['to'] in entrezToUniProt:
            entrezToUniProt[line['to']].append(line['from'])
        else:
            entrezToUniProt[line['to']] = [line['from']]

    for line in notMappedData:
        if line['from'] not in uniProtToEntrez:
            uniProtToEntrez[line['from']] = [line['to']]

    return uniProtToEntrez, entrezToUniProt


#get_uniprot_to_entrez_map()
jsonToDict()
