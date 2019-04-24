import os
from gensim.models import Word2Vec
from lib import node2vec
import datetime
import config
import numpy as np

'''
TODO: migrate to https://github.com/eliorc/node2vec
'''
def process(pathway_id, nx_G, args, gene_vec_conv=lambda x: x):
    '''Generates node2vec representation of genes for given pathway network
    Parameters
    ----------
    pathway_id: str
        pathway that is being processed, will be used to cache processed pathway
    nx_G: :networkx.classes.graph.Graph:
        pathway graph
    args: dict
        node2vec arguments see node2vec
    gene_vec_conv: function
        gene vector converter function
    '''
    # generate model
    for v1, v2 in nx_G.edges(): nx_G[v1][v2]['weight'] = 1
    G = node2vec.Graph(nx_G, is_directed=args.is_directed, p=args.p, q=args.q)
    G.preprocess_transition_probs()
    walks_sim = G.simulate_walks(num_walks=10, walk_length=80)
    walks = [list(map(str, walk)) for walk in walks_sim]
    # size=dimension, window=context_size, workers=num_of_parallel_workers, iter=num_epochs_in_sgd
    model = Word2Vec(walks, size=args.n2v_size, window=10, min_count=0, sg=1, workers=4, iter=1)

    # disabled backup because caller will handle it
    ## backup node2vec result
    # FNAME = '{}-p={:0.2f}-q={:0.2f}-dir={}-run={}-word2vec.csv'.format(pathway_id, args.p, args.q, args.is_directed, args.rid)
    # FPATH = os.path.join(config.data_dir, FNAME)
    # model.wv.save_word2vec_format(FPATH)

    # wv has:
    #   * index2entity: array of result index to node id map size N
    #   * vectors: array of result vectors size N x v2v_size
    # so convert them to map G (gene vector map)
    gene_vectors = {}
    for (eid, gene_vec) in zip(model.wv.index2entity, model.wv.vectors):
        gene_vectors[int(eid)] = gene_vec_conv(gene_vec)

    return gene_vectors
