import os
import config
import datetime
import requests


def initialize_node2vec():
    print('node2vec: checking')
    # load node2vec as python 3 module
    NODE2VEC_NAME = 'node2vec.py'
    NODE2VEC_PATH = os.path.join(config.lib_dir, NODE2VEC_NAME)
    if not os.path.exists(NODE2VEC_PATH):
        r = requests.get('https://raw.githubusercontent.com/aditya-grover/node2vec/master/src/node2vec.py')
        with open(NODE2VEC_PATH, 'w') as q:
            for line in r.text.split('\n'):
                line = line.replace('save_word2vec_format', 'wv.save_word2vec_format')
                if 'print' in line:
                    line = line.replace('print', 'if self.debug: print(') + ')'
                    # gensim version update
                if 'def __init__' in line:
                    line = f'{line[:-2]}, debug=False):\n\t\tself.debug = debug'
                q.write(f'{line}\n')
            q.write(f'\n# accessed:{datetime.datetime.now()}')
        print('node2vec: installed')
    else:
        print('node2vec: already exists skipping')

if not os.path.exists(config.lib_dir):
    print('libs_dir: creating')
    os.makedirs(config.lib_dir)

initialize_node2vec()

