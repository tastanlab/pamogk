#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import config
import requests
import os
import json
import networkx as nx
import pdb
from lib.sutils import *

HOST = 'http://www.ndexbio.org/v2'
NCI_USER_ID = '301a91c6-a37b-11e4-bda0-000c29202374'

DATA_ROOT = os.path.join(config.data_dir, 'cx')
safe_create_dir(DATA_ROOT)

def get_pathway_map():
    PATHWAY_LIST_PATH = os.path.join(DATA_ROOT, 'nci_pathway_list.json')

    if not os.path.exists(PATHWAY_LIST_PATH):
        try:
            url = '{}/user/{}/showcase'.format(HOST, NCI_USER_ID)
            log('Pathway map not found fetching from', url)
            r = requests.get(url)
            if r.status_code != 200:
                raise Exception('Failed to get pathway list:', r.reason)
        except requests.exceptions.ConnectionError as e:
            log('Could not connect to {} server. Please make sure you are connected to a network authorized to access {}'.format(HOST, HOST))
            raise Exception('Failed to connect to host={}'.format(HOST))

        with open(PATHWAY_LIST_PATH, 'w') as f:
            f.write(r.text)

    pathway_list = json.load(open(PATHWAY_LIST_PATH))

    pathway_map = {}
    for p in pathway_list:
        pathway_map[p['externalId']] = p
    return pathway_map

def _get_pathway_child(pathway_data, key):
    for d in pathway_data:
        if key in d:
            return d[key]
    return None

@timeit
def read_pathways():
    pathway_map = get_pathway_map()
    pw_map = {}
    pw_ids = pathway_map.keys()
    log('Pathway data_dir={}'.format(DATA_ROOT))
    for (ind, pw_id) in enumerate(pw_ids):
        log('Processing pathway {:3}/{}'.format(ind + 1, len(pw_ids)), end='\t')
        pw_data = read_single_pathway(pw_id, reading_all=True)
        pw_map[pw_id] = pw_data
    log()
    return pw_map

def read_single_pathway(pathway_id, reading_all=False):
    pend = '\r' if reading_all else '\n'
    pathway_map = get_pathway_map()
    if pathway_id not in pathway_map:
        raise Exception('Pathway not found in pathway list')

    PATHWAY_PATH = os.path.join(DATA_ROOT, pathway_id + '.cx')

    if not os.path.exists(PATHWAY_PATH):
        url = '{}/network/{}'.format(HOST, pathway_id)
        log('Pathway with pathway_id={} not found fetching from url={}'.format(pathway_id, url), end=pend, ts=not reading_all)
        r = requests.get(url)
        if r.status_code != 200:
            raise Exception('Failed to get pathway: {}-{}'.format(r.status_code, r.reason))

        with open(PATHWAY_PATH, 'w') as f:
            f.write(r.text)
    else:
        log('Pathway with pathway_id={} retrieved from local data dir'.format(pathway_id), end=pend, ts=not reading_all)

    pathway_data = json.load(open(PATHWAY_PATH))

    G = nx.Graph() # initialize empty graph

    # get node map
    node_list = _get_pathway_child(pathway_data, 'nodes')
    nodes = {}
    for n in node_list: nodes[n['@id']] = n

    # load node attributes
    node_attributes = _get_pathway_child(pathway_data, 'nodeAttributes')
    attr_dict = {}
    for attr in node_attributes:
        nid = attr['po']
        if nid not in attr_dict:
            attr_dict[nid] = {}
        attr_dict[nid][attr['n']] = attr['v']

    # load cartesian coordinates
    coords = _get_pathway_child(pathway_data, 'cartesianLayout')
    for coord in coords:
        nid = coord['node']
        if nid not in attr_dict:
            attr_dict[nid] = {}
        attr_dict[nid]['x'] = coord['x']
        attr_dict[nid]['y'] = coord['y']

    # add nodes to graph
    # NOTE networkx graphs only allow alphanumeric characters as attribute names no - or _
    for nid in nodes:
        n = nodes[nid]
        attrs = attr_dict[nid]
        attrs['n'] = n['n']
        if 'r' in n:
            if 'alias' not in attrs: attrs['alias'] = [n['r']]
            else: attrs['alias'].append(n['r'])
        if 'alias' in attrs: # create attribute for ids only
            tmp = [a.split(':') for a in attrs['alias']]
            attrs['uniprot-ids'] = [s[1] for s in tmp if len(s) > 1 and len(s[1]) > 0]
        G.add_node(nid, **attrs)

    # get edge map
    edge_list = _get_pathway_child(pathway_data, 'edges')
    edges = {}
    for e in edge_list:
        edges[e['@id']] = e

    # load edge attributes
    edge_attributes = _get_pathway_child(pathway_data, 'edgeAttributes')
    attr_dict = {}
    for attr in edge_attributes:
        eid = attr['po']
        if eid not in attr_dict:
            attr_dict[eid] = {}
        attr_dict[eid][attr['n']] = attr['v']
    # add edges to graph
    for eid in edges:
        e = edges[eid]
        attrs = attr_dict[eid]
        attrs['i'] = e['i']
        G.add_edge(e['s'], e['t'], **attrs)

    return G

if __name__ == '__main__':
    G = read_pathways()
    pdb.set_trace()
