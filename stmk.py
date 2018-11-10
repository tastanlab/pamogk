#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Usage: python3 stmk.py hsa04151

Packages required:
requests
biopython
networkx
numpy
gensim
plotly
'''

import requests
import sys
import os
import math
import networkx as nx
import numpy as np
import config
import argparse
# from pathway_reader import kgml_converter
from pathway_reader import cx_pathway_reader
from pathway_reader import network_plotter
from structural_processor import node2vec_processor
from synthetic_experiments import cell_survival_group_kegg
from lib import tsne
from sklearn.svm import SVC
from gene_mapper import uniprot_mapper

sys.path.append(config.root_dir)

# import plotly.offline as py
import plotly.plotly as py
import plotly.graph_objs as go

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('pathways', metavar='pathway-id', type=str, nargs='+', help='pathway ID list', default=['hsa04151'])
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float, help='Node2Vec feature space size', default=128)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed', default=False)
parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments', default=1000)
parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float, help='Surviving patient percentage in range [0, 1]', default=0.9)
parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]', default=0.4)

args = parser.parse_args()
print('Running args:', args)

# To fallback to python debug console
# import pdb; pdb.set_trace()

# get pathway id from arguments if given
pathway_id = args.pathways[0]
pathway_id = '8bbf39aa-6193-11e5-8ac5-06603eb7f303'

OUT_FILENAME = os.path.join(config.data_dir, '{}-p={:0.2f}-q={:0.2f}-undirected-run={}'.format(pathway_id, args.p, args.q, args.rid))

# read kgml in format of network
# nx_G, entries, relations = kgml_converter.KGML_to_networkx_graph(pathway_id, is_directed=args.is_directed)
# all_pws = cx_pathway_reader.read_pathways()
uniprot_mapper.get_uniprot_to_entrez_map()
# import pdb; pdb.set_trace()

sys.exit(0)
nx_G = cx_pathway_reader.read_single_pathway(pathway_id)
network_plotter.plot(nx_G, title='Pathway Graph for {}'.format(pathway_id))
# network_plotter.plot(nx_G, title='TSNE Representation of {} p={:0.2f} q={:0.2f} run={}'.format(pathway_id, args.p, args.q, args.rid))

# run node2vec to get feature representations
gene_vectors = node2vec_processor.process(pathway_id, nx_G, args)
hnames = np.array([nx_G.node[int(eid)]['hname'] for eid in gene_vectors])

'''
patients = cell_survival_group_kegg.generate_patients(G=nx_G, num_pat=args.num_pat, surv_dist=args.surv_dist, mut_dist=args.mut_dist)

center_product_kernel.calculate_S_and_P(patients, gene_vectors)
center_product_kernel.test_accr(patients)

sys.exit(0)
'''

# visualize node2vec nodes
# run tsne
Y = tsne.tsne(X=gene_vectors.values(), no_dims=2, initial_dims=50, perplexity=30.0)
# np.savetxt(pathway_id + '-tsne-hnames.csv', hnames)
TSNE_OUT_PATH = OUT_FILENAME + '-tsne.csv'
np.savetxt(TSNE_OUT_PATH, X)

# cluster using dbscan
from sklearn import cluster

# clt = cluster.SpectralClustering(n_clusters=3, eigen_solver='arpack', affinity="nearest_neighbors").fit(X)
clt = cluster.AffinityPropagation().fit(X)
# clt = cluster.DBSCAN(eps=0.3, min_samples=10).fit(X)
clabels = clt.labels_
# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(clabels)) - (1 if -1 in clabels else 0)
print('n_clusters =', n_clusters_)

# import colorlover as cl
# bupu = cl.scales['9']['seq']['BuPu']
# palette = cl.interp( bupu, max(n_clusters_, 10))
# colors = ['rgb(0, 0, 0)' if k == -1 else palette[k] for k in labels]

labels = ['{}-{}'.format(t.split(',')[0], gid) for t, gid in zip(hnames, clabels)]
labels = ['?' + l if l[0] == '-' else l for l in labels]
print(labels)

fig = go.Figure(
    data = [
        go.Scatter(
        x=Y[:, 0],
        y=Y[:, 1],
        mode='markers+text',
        text=labels,
        textposition='bottom center',
        marker=go.scatter.Marker(color=clabels, colorscale='Viridis', showscale=True))
    ],
    layout = go.Layout(
        title= 'TSNE Representation of {} p={:0.2f} q={:0.2f} run={}'.format(pathway_id, args.p, args.q, args.rid),
        hovermode= 'closest',
        xaxis= dict(
            ticklen= 5,
            zeroline= False,
            gridwidth= 2,
        ),
        yaxis = dict(
            ticklen=5,
            gridwidth=2,
        ),
        showlegend=False,
        width=1200,
        height=900,
    ))
import plotly.offline as pyoff
pyoff.plot(fig, filename=OUT_FILENAME + '-plot.html', auto_open=False)
py.plot(fig, filename=OUT_FILENAME + '-plot.html', auto_open=False)

import traceback
for i in range(5):
    try:
        IMG_PATH = OUT_FILENAME + '-plot.png'
        py.image.save_as(fig, filename=IMG_PATH)
        print('Saved plot image to:', IMG_PATH)
        break
    except Exception as e:
        traceback.print_exc(e)
