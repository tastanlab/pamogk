import datetime
from pathlib import Path

import requests

from .sutils import safe_create_dir
from .. import config


def initialize_node2vec():
    LIB_DIR = Path(__file__).resolve().parent / 'lib'

    if not LIB_DIR.exists():
        print('libs_dir: creating')
        safe_create_dir(LIB_DIR)

    print('node2vec: checking')
    # load node2vec as python 3 module
    node2_vec_name = 'node2vec.py'
    node2_vec_path = config.LIB_DIR / node2_vec_name
    if not node2_vec_path.exists():
        r = requests.get('https://raw.githubusercontent.com/aditya-grover/node2vec/master/src/node2vec.py')
        with open(node2_vec_path, 'w') as q:
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


initialize_node2vec()
