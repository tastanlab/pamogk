#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
import and use methods

Packages required:
requests
biopython
networkx
numpy
gensim
plotly
'''
import os
import requests
import config
import networkx as nx
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.KGML.KGML_pathway import Relation
from . import kgml_pathway_reader

def KGML_to_networkx_graph(pathway_id, is_directed, entries=None, relations=None):
    '''
    @param is_directed forced to prevent parsing the graph wrongly
    '''
    print('Converting KGML pathway to networkx pathway:', pathway_id)
    if entries is None or relations is None:
        entries, relations = kgml_pathway_reader.get_pathway_KGML(pathway_id)

    # validate data dir
    if not os.path.exists(config.data_dir): os.makedirs(config.data_dir)
    # convert kgml data to networkx graph
    nodes = [(eid, {'name': entries[eid].name, 'hname': entries[eid].graphics[0].name, 'eid': eid }) for eid in entries]
    edges = [(r.entry1.id, r.entry2.id, {'weight': 1}) for r in relations]
    nx_G = nx.Graph()
    nx_G.add_nodes_from(nodes)
    nx_G.add_edges_from(edges)
    # remove 0 degree info nodes
    # nx_G.remove_nodes_from([n for n in nx.isolates(nx_G)])
    print('Nodes with no edges:', len([n for n in nx.isolates(nx_G)]))
    path = os.path.join(config.data_dir, pathway_id + '.gml')
    print('Saving GML to path:', path)
    nx.write_gml(nx_G, path)

    print('Finished conversion to networkx graph pathway:', pathway_id)
    if is_directed:
        return nx_G, entries, relations
    else:
        return nx_G.to_undirected(), entries, relations