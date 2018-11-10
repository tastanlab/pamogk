#!/usr/bin/env python3
import requests
import config
import numpy as np
from pathway_reader import cx_pathway_reader

def get_all_nodes(all_pws_map):
    return np.hstack([[g.nodes[n] for n in g.nodes()] for g in all_pws_map.values()])

def get_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.hstack([n['alias'] for n in all_nodes if 'alias' in n])

def get_non_alias_set(all_pws_map):
    all_nodes = get_all_nodes(all_pws_map)
    return np.hstack([n['n'] for n in all_nodes if 'alias' not in n])

def fetch_uniprot_mapping(alias_list, fr='ACC,Id', to='P_ENTREZGENEID'):
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
            'uploadQuery': alias_list[0],
            'url2': '',
            'from': fr,
            'to': to,
            'taxon': '',
        },
        files = { 'file': None })
    MAPPING_ID = None
    if len(r.history) > 0: MAPPING_ID = r.history[0].headers['location'].split('/')[-1]
    if MAPPING_ID is None:
        raise Exception('Could not find mapping for aliases')
    r1 = requests.get('https://www.uniprot.org/mapping/%s.tab' % MAPPING_ID)
    print(r1.text)
    r2 = requests.get('https://www.uniprot.org/mapping/%s.not' % MAPPING_ID)
    print(r2.text)
    import pdb; pdb.set_trace()

def fetch_uniprot_to_entrez(alias_list):
    uniprot_alias_list = [a[1] for a in [a.split('uniprot knowledgebase:')for a in alias_list] if len(a) > 1]
    return fetch_uniprot_mapping(uniprot_alias_list)

def get_uniprot_to_entrez_map(all_pws_map = None):
    if all_pws_map is None:
        all_pws_map = cx_pathway_reader.read_pathways()
    all_alias = get_alias_set(all_pws_map)
    data = fetch_uniprot_to_entrez(all_alias)
    path = os.path.join(config.data_dir, 'gene-map-uniprot-entrez.json')
    with open(path) as f:
        json.dump(data, f)
