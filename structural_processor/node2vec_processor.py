import os
import shutil
from subprocess import call

def initailize_node2vec():
    # load node2vec as python 3 module
    NODE2VEC_NAME = 'node2vec.py'
    NODE2VEC_PATH = os.path.join(os.path.realpath(os.path.dirname(__file__)), NODE2VEC_NAME)
    NODE2VEC_GIT_DIR = 'node2vec'
    if not os.path.exists(NODE2VEC_PATH):
        call(['git', 'clone', 'https://github.com/aditya-grover/node2vec', NODE2VEC_GIT_DIR])
        try:
            with open('node2vec/src/' + NODE2VEC_NAME) as f, open(NODE2VEC_PATH, 'w') as q:
                for line in f.readlines():
                    if 'print' in line:
                        line = line.replace('print', 'print(')[:-1] + ')\n'
                    q.write(line)
            shutil.rmtree(NODE2VEC_GIT_DIR)
        except Exception as e:
            os.remove(NODE2VEC_PATH)
            raise e

    from . import node2vec

def process(args):
    initailize_node2vec()
    # generate model
    G = node2vec.Graph(nx_G, is_directed=args.is_directed, p=args.p, q=args.q)
    G.preprocess_transition_probs()
    walks_sim = G.simulate_walks(num_walks=10, walk_length=80)
    walks = [list(map(str, walk)) for walk in walks_sim]
    # size=dimension, window=context_size, workers=num_of_parallel_workers, iter=num_epochs_in_sgd
    model = Word2Vec(walks, size=128, window=10, min_count=0, sg=1, workers=4, iter=1)
    WORD2VEC_PATH = OUT_FILENAME + '-word2vec.csv'
    model.wv.save_word2vec_format(WORD2VEC_PATH)
    return

