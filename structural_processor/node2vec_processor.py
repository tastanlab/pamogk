import os
from gensim.models import Word2Vec
from lib import node2vec
import datetime
import config
from pathway_reader import cx_pathway_reader
from structural_processor import node2vec_processor
import argparse

def process(pathway_id, nx_G, args):
    '''
    Generates node2vec representation of genes for given pathway network
    '''
    # generate model
    print('Calculating node2vec with params:', args)
    for v1, v2 in nx_G.edges(): nx_G[v1][v2]['weight'] = 1
    G = node2vec.Graph(nx_G, is_directed=args.is_directed, p=args.p, q=args.q)
    G.preprocess_transition_probs()
    walks_sim = G.simulate_walks(num_walks=10, walk_length=80)
    walks = [list(map(str, walk)) for walk in walks_sim]
    # size=dimension, window=context_size, workers=num_of_parallel_workers, iter=num_epochs_in_sgd
    model = Word2Vec(walks, size=args.n2v_size, window=10, min_count=0, sg=1, workers=4, iter=1)
    # backup node2vec result
    FNAME = '{}-p={:0.2f}-q={:0.2f}-dir={}-run={}-word2vec.csv'.format(pathway_id, args.p, args.q, args.is_directed, args.rid)
    FPATH = os.path.join(config.data_dir+"\\node2vec", FNAME)
    model.wv.save_word2vec_format(FPATH)
    # wv has:
    #   * index2entity: array of result index to node id map size N
    #   * vectors: array of result vectors size N x v2v_size
    # so convert them to map G (gene vector map)
    gene_vectors = {}
    for (eid, gene_vec) in zip(model.wv.index2entity, model.wv.vectors):
        gene_vectors[int(eid)] = gene_vec

    return gene_vectors

##To create vectors for all pathways
def processAllPathways():
    parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
    parser.add_argument('--pathways', metavar='pathway-id', type=str, nargs='+', help='pathway ID list',
                        default=['hsa04151'])
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


    pathways = cx_pathway_reader.get_pathway_map()
    for aPathway in pathways:
        nx_G = cx_pathway_reader.read_single_pathway(aPathway)
        gene_vec_map = node2vec_processor.process(aPathway, nx_G, args)

