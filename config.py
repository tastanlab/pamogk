from lib.sutils import *

root_dir = os.path.realpath(os.path.dirname(__file__))
lib_dir = os.path.join(root_dir, 'lib')
data_dir = os.path.join(root_dir, 'data')

safe_create_dir(lib_dir)
safe_create_dir(data_dir)
