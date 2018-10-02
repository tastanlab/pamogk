import os
from gensim.models import Word2Vec
from lib import node2vec
import datetime
import config

def process(pathway_id, nx_G, args):
    # generate model
    G = node2vec.Graph(nx_G, is_directed=args.is_directed, p=args.p, q=args.q)
    G.preprocess_transition_probs()
    walks_sim = G.simulate_walks(num_walks=10, walk_length=80)
    walks = [list(map(str, walk)) for walk in walks_sim]
    # size=dimension, window=context_size, workers=num_of_parallel_workers, iter=num_epochs_in_sgd
    model = Word2Vec(walks, size=128, window=10, min_count=0, sg=1, workers=4, iter=1)
    OUT_FILENAME = os.path.join(config.data_dir, '{}-p={:0.2f}-q={:0.2f}-undirected-run={}'.format(pathway_id, args.p, args.q, args.rid))
    WORD2VEC_PATH = OUT_FILENAME + '-word2vec.csv'
    model.wv.save_word2vec_format(WORD2VEC_PATH)
    print(model.wv)
    return

