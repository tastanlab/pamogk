#!/usr/bin/env python3
import json
import pdb

import numpy as np
import requests

from .. import config
from ..pathway_reader import cx_pathway_reader

UNIPROT_DOMAIN = 'www.uniprot.org'
UNIPROT_HOST = f'https://{UNIPROT_DOMAIN}'


def get_all_nodes(all_pws_map):
    return np.hstack([[g.nodes[n] for n in g.nodes()] for g in all_pws_map.values()])


def get_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.unique(np.hstack([n['alias'] for n in all_nodes if 'alias' in n]))


def get_non_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.unique(np.hstack([n['n'] for n in all_nodes if 'alias' not in n]))


def fetch_uniprot_mapping(alias_list, fr='ACC,ID', to='P_ENTREZGENEID'):
    print('Requesting mapping for', len(alias_list), 'genes')
    alias_list_str = ",".join(alias_list)
    path = config.DATA_DIR / 'all_aliases.txt'
    with open(path, 'w') as f:
        f.write(alias_list_str)
    url = f'{UNIPROT_HOST}/uploadlists/'
    r = requests.post(url,
                      headers={
                          'authority': UNIPROT_DOMAIN,
                          'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36',
                          'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
                          'upgrade-insecure-requests': '1',
                          'referer': url,
                      },
                      data={
                          'landingPage': 'false',
                          'jobId': '',
                          'uploadQuery': alias_list_str,
                          'url2': '',
                          'from': fr,
                          'to': to,
                          'taxon': '',
                      })
    mapping_id = None
    if len(r.history) > 0:
        mapping_id = r.history[0].headers['location'].split('/')[-1]

    if r.status_code != requests.codes.ok:
        raise Exception(f'Failed to request mapping data with status={r.status_code} and reason={r.text}')

    if mapping_id is None:
        raise Exception('Could not find mapping for aliases')

    print('Received mapping with id:', mapping_id)

    r = requests.get('%s/mapping/%s.tab' % (UNIPROT_HOST, mapping_id))
    if r.status_code != requests.codes.ok:
        raise Exception(f'Failed to get mapping data with status={r.status_code} and reason={r.text}')

    rows = [r.split('\t') for r in r.text.split('\n')[1:-1]]
    mapped = [{'from': r[0], 'to': r[1]} for r in rows]

    r = requests.get('%s/mapping/%s.not' % (UNIPROT_HOST, mapping_id))
    if r.status_code != requests.codes.ok:
        raise Exception(f'Failed to get missing data with status={r.status_code} and reason={r.text}')

    missing = [r for r in r.text.split('\n')[1:-1]]

    return mapped, missing


def fetch_uniprot_to_entrez(alias_list):
    uniprot_alias_list = [a[1] for a in [a.split('uniprot knowledgebase:') for a in alias_list] if len(a) > 1]
    return fetch_uniprot_mapping(uniprot_alias_list)


def get_uniprot_to_entrez_map(all_pws_map=None):
    if all_pws_map is None:
        all_pws_map = cx_pathway_reader.read_pathways()

    all_alias = get_alias_set(all_pws_map)

    mapped, missing = fetch_uniprot_to_entrez(all_alias)

    json_data = {
        'mapped': mapped,
        'missing': missing,
    }

    path = config.DATA_DIR / 'gene-map-uniprot-entrez.json'
    with open(path, 'w') as f:
        json.dump(json_data, f)


def json_to_dict():
    path = config.DATA_DIR / 'gene-map-uniprot-entrez.json'
    with open(path, 'r') as f:
        data = json.load(f)

    mapped_data = data['mapped']
    missing_data = data['missing']

    uni_prot_to_entrez = {}
    entrez_to_uni_prot = {}

    for line in mapped_data:
        safe_list_value_append(uni_prot_to_entrez, line['from'], line['to'])
        safe_list_value_append(entrez_to_uni_prot, line['to'], line['from'])

    for gene in missing_data:
        if gene not in uni_prot_to_entrez:
            uni_prot_to_entrez[gene] = []

    return uni_prot_to_entrez, entrez_to_uni_prot


def safe_list_value_append(d, k, v):
    if k not in d:
        d[k] = [v]
    else:
        d[k].append(v)


if __name__ == '__main__':
    get_uniprot_to_entrez_map()
    json_to_dict()
    pdb.set_trace()
