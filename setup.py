import os
import shutil
import config
import datetime
import requests
from subprocess import call
from zipfile import ZipFile
from io import BytesIO
from distutils.core import setup
import setuptools

def initialize_node2vec():
    print('node2vec: checking')
    # load node2vec as python 3 module
    NODE2VEC_NAME = 'node2vec.py'
    NODE2VEC_PATH = os.path.join(config.lib_dir, NODE2VEC_NAME)
    NODE2VEC_GIT_DIR = 'node2vec'
    if not os.path.exists(NODE2VEC_PATH):
        call(['git', 'clone', 'https://github.com/aditya-grover/node2vec', NODE2VEC_GIT_DIR])
        try:
            with open('node2vec/src/' + NODE2VEC_NAME) as f, open(NODE2VEC_PATH, 'w') as q:
                for line in f.readlines():
                    if 'print' in line:
                        line = line.replace('print', 'print(')[:-1] + ')\n'
                    q.write(line)
                q.write('# accessed:' + str(datetime.datetime.now()))
            shutil.rmtree(NODE2VEC_GIT_DIR)
            print('node2vec: installed')
        except Exception as e:
            os.remove(NODE2VEC_PATH)
            raise e
    else:
        print('node2vec: already exists skipping')

def initialize_tsne():
    print('tsne: checking')
    # load tsne
    TSNE_NAME = 'tsne.py'
    TSNE_PATH = os.path.join(config.lib_dir, TSNE_NAME)
    if not os.path.exists(TSNE_PATH):
        r = requests.get('https://lvdmaaten.github.io/tsne/code/tsne_python.zip')
        if r.status_code == 200:
            zipfile = ZipFile(BytesIO(r.content))
            print('Files in zip:\n\t', '\n\t'.join([fn for fn in zipfile.namelist() if fn[:8] != '__MACOSX']), sep='')
            try:
                with open(TSNE_PATH, 'w') as f:
                    for line in zipfile.open('tsne_python/' + TSNE_NAME).readlines():
                        lstr = line.decode('utf-8')
                        # remove demo parts
                        if '__main__' in lstr:
                            break
                        if 'pylab' not in lstr:
                            f.write(lstr)
                    f.write('# accessed:' + str(datetime.datetime.now()))
                print('tsne: installed')
            except Exception as e:
                os.remove(TSNE_PATH)
                raise e
    else:
        print('tsne: already exists skipping')

if not os.path.exists(config.lib_dir):
    print('libs_dir: creating')
    os.makedirs(config.lib_dir)

initialize_node2vec()
initialize_tsne()

# load requirements from file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.txt') as f:
    long_description = f.read()

setup(
    name='smSPK',
    version='0.1.0',
    packages=['smspk',],
    install_requires=requirements,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=long_description,
)
