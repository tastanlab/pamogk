import time
import os
import sys
import pdb
import numpy as np

def timeit(f):
    '''Timing decorator for functions. Just add @timeit to start and function
    will be timed. It will write starting and ending times

    Parameters
    ----------
    f : function
        decorator takes the function as parameter

    Returns
    -------
    mixed
        return value of function itself

    Raises
    ------
    Error
        when any error thrown by called function does not catch it
    '''

    def wrapper(*args, **kwargs):
        log("Started:", f.__qualname__)
        t = time.time()
        res = f(*args, **kwargs)
        log("Finished: {} elapsed: {:.2f}s".format(f.__qualname__, time.time() - t))
        return res
    return wrapper

def safe_create_dir(d):
    if not os.path.exists(d):
        log('Dir not found creating:', d)
        os.makedirs(d)

log_f = None
log_p = None

def change_log_path(path):
    global log_p, log_f
    if log_p == path: return
    if log_f: log_f.close()
    log_p = path
    d = os.path.basename(path)
    safe_create_dir(d)
    log_f = open(path, 'a')
    log('Initialized log_path:', path)

def log(*args, **kwargs):
    print(*args, **kwargs)
    if log_f:
        print(*args, **kwargs, file=log_f)
        log_f.flush()

def exists_np_data(path):
    return os.path.exists(path + '.npz')

def save_np_data(path, *args, **kwargs):
    if len(kwargs) == 0: kwargs = { 'data': args[0] }
    np.savez_compressed(path + '.npz', **kwargs)

def load_np_data(path, key='data'):
    return np.load(path + '.npz')['data']