#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: python3 stmk.py hsa04151

See readme
"""

import argparse

import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
from sklearn import cluster
from sklearn.manifold import TSNE

from pamogk import config
from pamogk.data_processor import node2vec_processor
from pamogk.lib.sutils import print_args
from pamogk.pathway_reader import cx_pathway_reader
from visualizations import network_plotter

parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('pathways', metavar='pathway-id', type=str, nargs='+', help='pathway ID list',
                    default=['8bbf39aa-6193-11e5-8ac5-06603eb7f303'])
parser.add_argument('--debug', action='store_true', dest='debug', help='Enable Debug Mode')
parser.add_argument('--node2vec-p', '-p', metavar='p', dest='p', type=float, help='Node2Vec p value', default=1)
parser.add_argument('--node2vec-q', '-q', metavar='q', dest='q', type=float, help='Node2Vec q value', default=1)
parser.add_argument('--node2vec-size', '-n', metavar='node2vec-size', dest='n2v_size', type=float,
                    help='Node2Vec feature space size', default=128)
parser.add_argument('--run-id', '-r', metavar='run-id', dest='rid', type=str, help='Run ID', default=None)
parser.add_argument('--directed', '-d', dest='is_directed', action='store_true', help='Is graph directed',
                    default=False)
parser.add_argument('--num-pat', dest='num_pat', type=int, help='Number of Patients for Synthetic Experiments',
                    default=1000)
parser.add_argument('--surv-dist', '-s', dest='surv_dist', type=float,
                    help='Surviving patient percentage in range [0, 1]', default=0.9)
parser.add_argument('--mut-dist', '-m', dest='mut_dist', type=float, help='Mutated gene percentage in range [0, 1]',
                    default=0.4)

args = parser.parse_args()
print_args(args)

# To fallback to python debug console
# import pdb; pdb.set_trace()

# get pathway id from arguments if given
pathway_id = args.pathways[0]

OUT_FILENAME = config.DATA_DIR / f'{pathway_id}-p={args.p:0.2f}-q={args.q:0.2f}-undirected-run={args.rid}'

# read kgml in format of network
# nx_g, entries, relations = kgml_converter.KGML_to_networkx_graph(pathway_id, is_directed=args.is_directed)
# all_pws = cx_pathway_reader.read_pathways()
# uniprot_mapper.get_uniprot_to_entrez_map()

# sys.exit(0)
nx_G = cx_pathway_reader.read_single_pathway(pathway_id)
network_plotter.plot(nx_G, title=f'Pathway Graph for {pathway_id}')
# network_plotter.plot(nx_G, title=f'TSNE Representation of {pathway_id} p={args.p:0.2f} q={args.q:0.2f} run={args.rid}')

# run node2vec to get feature representations
gene_vec_map = node2vec_processor.process(nx_G, args)
hnames = np.array([nx_G.nodes[int(eid)]['n'] for eid in gene_vec_map])

'''
patients = cell_survival_group_kegg.generate_patients(G=nx_g, num_pat=args.num_pat, surv_dist=args.surv_dist, mut_dist=args.mut_dist)

center_product_kernel.calculate_S_and_P(patients, gene_vec_map)
center_product_kernel.test_accr(patients)

sys.exit(0)
'''

# visualize node2vec nodes
# run tsne
gene_vectors = np.array([x for x in gene_vec_map.values()])
T = TSNE(n_components=2, perplexity=30.0).fit_transform(gene_vectors)
# np.savetxt(pathway_id + '-tsne-hnames.csv', hnames)
TSNE_OUT_PATH = OUT_FILENAME + '-tsne.csv'
np.savetxt(TSNE_OUT_PATH, T)

# cluster using dbscan

# clt = cluster.SpectralClustering(n_clusters=3, eigen_solver='arpack', affinity='nearest_neighbors').fit(gene_vectors)
clt = cluster.AffinityPropagation().fit(gene_vectors)
# clt = cluster.DBSCAN(eps=0.3, min_samples=10).fit(gene_vectors)
clabels = clt.labels_
# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(clabels)) - (1 if -1 in clabels else 0)
print('n_clusters =', n_clusters_)

# import colorlover as cl
# bupu = cl.scales['9']['seq']['BuPu']
# palette = cl.interp( bupu, max(n_clusters_, 10))
# colors = ['rgb(0, 0, 0)' if k == -1 else palette[k] for k in labels]

labels = [f'{t.split(",")[0]}-{gid}' for t, gid in zip(hnames, clabels)]
labels = [f'?{line}' if line[0] == '-' else line for line in labels]
print(labels)

fig = go.Figure(
    data=[
        go.Scatter(
            x=T[:, 0],
            y=T[:, 1],
            mode='markers+text',
            text=labels,
            textposition='bottom center',
            marker=go.scatter.Marker(color=clabels, colorscale='Viridis', showscale=True))
    ],
    layout=go.Layout(
        title=f'TSNE Representation of {pathway_id} p={args.p:0.2f} q={args.q:0.2f} run={args.rid}',
        hovermode='closest',
        xaxis=dict(
            ticklen=5,
            zeroline=False,
            gridwidth=2,
        ),
        yaxis=dict(
            ticklen=5,
            gridwidth=2,
        ),
        showlegend=False,
        width=1200,
        height=900,
    ))
pio.write_html(fig, f'{OUT_FILENAME}-plot.html',
               include_plotlyjs='directory',
               auto_open=True)
'''
# save library to plotly cloud for online sharing
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
'''
