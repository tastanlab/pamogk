import os

root_dir = os.path.realpath(os.path.dirname(__file__))
lib_dir = os.path.join(root_dir, 'lib')
data_dir = os.path.join(root_dir, 'data')

def create_if_not_exists(dirpath):
    if not os.path.exists(dirpath):
        print('Initializing dir:', dirpath)
        os.makedirs(dirpath)

create_if_not_exists(lib_dir)
create_if_not_exists(data_dir)
