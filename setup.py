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

if not os.path.exists(config.lib_dir):
    print('libs_dir: creating')
    os.makedirs(config.lib_dir)

initialize_node2vec()

# load requirements from file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md') as f:
    long_description = f.read()

setup(
    name='pamogk',
    version='0.1.0',
    packages=['pamogk',],
    install_requires=requirements,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=long_description,
)
