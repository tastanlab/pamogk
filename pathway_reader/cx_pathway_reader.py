#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import config
import requests
import os
import json
import networkx as nx

HOST = 'http://www.ndexbio.org/v2'
NCI_USER_ID = '301a91c6-a37b-11e4-bda0-000c29202374'

def get_pathway_map():
    PATHWAY_LIST_PATH = os.path.join(config.data_dir, 'nci_pathway_list.json')

    if not os.path.exists(PATHWAY_LIST_PATH):
        r = requests.get('{}/user/{}/showcase'.format(HOST, NCI_USER_ID))
        if r.status_code != 200:
            raise Exception('Failed to get pathway list:', r.reason)

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

def read_pathways(pathway_id):
    pathway_map = get_pathway_map()
    if pathway_id not in pathway_map:
        raise Exception('Pathway not found in pathway list')

    PATHWAY_PATH = os.path.join(config.data_dir, pathway_id + '.cx')

    if not os.path.exists(PATHWAY_PATH):
        r = requests.get('{}/network/{}'.format(HOST, pathway_id))
        if r.status_code != 200:
            raise Exception('Failed to get pathway: {}-{}'.format(r.status_code, r.reason))

        with open(PATHWAY_PATH, 'w') as f:
            f.write(r.text)

    pathway_data = json.load(open(PATHWAY_PATH))

    G = nx.Graph()
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
        nodes[nid]['x'] = coord['x']
        nodes[nid]['y'] = coord['y']

    # add nodes to graph
    for nid in nodes:
        n = nodes[nid]
        G.add_node(nid, attr_dict=attr_dict[nid])

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

    import pdb; pdb.set_trace()
    return G
