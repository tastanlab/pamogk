#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Usage: python3 read-kgml.py hsa04151

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
from io import BytesIO
from zipfile import ZipFile
import networkx as nx
from gensim.models import Word2Vec
import numpy as np
from pathway_reader import kgml_converter
from structural_processor import node2vec_processor
import config
import argparse

# import plotly.offline as py
import plotly.plotly as py
import plotly.graph_objs as go

parser = argparse.ArgumentParser(description='Insert keywords from csv to account')
parser.add_argument('pathways', metavar='pathway-id', type=str, nargs='+', help='pathway ID list', default=['hsa04151'])
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed', default=False)

args = parser.parse_args()
print(args)

# get pathway id from arguments if given
pathway_id = args.pathways[0]

OUT_FILENAME = os.path.join(config.data_dir, '{}-p={:0.2f}-q={:0.2f}-undirected-run={}'.format(pathway_id, args.p, args.q, args.rid))

nx_G = kgml_converter.KGML_to_networkx_graph(pathway_id, is_directed=args.is_directed)

node2vec_features = node2vec_processor.process(nx_G)

# load tsne
TSNE_PATH = 'tsne.py'
if not os.path.exists(TSNE_PATH):
    r = requests.get('https://lvdmaaten.github.io/tsne/code/tsne_python.zip')
    if r.status_code == 200:
        zipfile = ZipFile(BytesIO(r.content))
        print('Files in zip:\n\t', '\n\t'.join([fn for fn in zipfile.namelist() if fn[:8] != '__MACOSX']), sep='')
        try:
            with open(TSNE_PATH, 'w') as f:
                for line in zipfile.open('tsne_python/' + TSNE_PATH).readlines():
                    lstr = line.decode('utf-8')
                    # remove demo parts
                    if '__main__' in lstr:
                        break
                    if 'pylab' not in lstr:
                        f.write(lstr)
        except Exception as e:
            os.remove(TSNE_PATH)
            raise e

import tsne

# visualize node2vec nodes
X = [list(map(float, line.strip().split(' '))) for line in open(WORD2VEC_PATH).readlines()[1:]]
hnames = np.array([nx_G.node[r[0]]['hname'] for r in X])
X = np.array([r[1:] for r in X])
# run tsne
Y = tsne.tsne(X=X, no_dims=2, initial_dims=50, perplexity=30.0)
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
