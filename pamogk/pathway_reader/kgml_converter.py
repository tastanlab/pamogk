#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import networkx as nx

from . import kgml_pathway_reader
from .. import config
from ..lib.sutils import *


@timeit
def KGML_to_networkx_graph(pathway_id, is_directed, entries=None, relations=None):
    """
    @param is_directed forced to prevent parsing the graph wrongly
    """
    print('Converting KGML pathway to networkx pathway:', pathway_id)
    if entries is None or relations is None:
        entries, relations = kgml_pathway_reader.get_pathway_kgml(pathway_id)

    # validate data dir
    safe_create_dir(config.DATA_DIR)
    # convert kgml data to networkx graph
    nodes = [(eid, {'name': entries[eid].name, 'hname': entries[eid].graphics[0].name, 'eid': eid }) for eid in entries]
    edges = [(r.entry1.id, r.entry2.id, {'weight': 1}) for r in relations]
    nx_G = nx.Graph()
    nx_G.add_nodes_from(nodes)
    nx_G.add_edges_from(edges)
    # remove 0 degree info nodes
    # nx_g.remove_nodes_from([n for n in nx.isolates(nx_g)])
    print('Nodes with no edges:', len([n for n in nx.isolates(nx_G)]))
    path = config.DATA_DIR / f'{pathway_id}.gml'
    print('Saving GML to path:', path)
    nx.write_gml(nx_G, path)

    print('Finished conversion to networkx graph pathway:', pathway_id)
    if is_directed:
        return nx_G, entries, relations
    else:
        return nx_G.to_undirected(), entries, relations
