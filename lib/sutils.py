import csv
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ISO_FORMAT = '%Y-%m-%dT%H:%M:%S.%f'


def timeit(f):
    """Timing decorator for functions. Just add @timeit to start and function
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
    """

    def wrapper(*args, **kwargs):
        log('Started:', f.__qualname__)
        t = time.time()
        res = f(*args, **kwargs)
        log(f'Finished: {f.__qualname__} elapsed: {time.time() - t:.2f}s')
        return res

    return wrapper


def safe_create_dir(d: Path):
    """
    Uses new pathlib
    Parameters
    ----------
    d: :obj:`pathlib.Path`
    """
    if not d.exists():
        log('Dir not found creating:', d)
        d.mkdir(parents=True)


def ensure_file_dir(file_path: Path):
    """
    Uses new pathlib
    Parameters
    ----------
    file_path: :obj:`pathlib.Path`
    """
    safe_create_dir(file_path.parent)


log_f = None
log_p = None


def change_log_path(path):
    global log_p, log_f
    if log_p == path:
        return
    if log_f:
        log_f.close()
    log_p = path
    ensure_file_dir(path)
    log_f = open(path, 'a')
    log('Initialized log_path:', path)


def logr(*args, **kwargs):
    log(*args, **kwargs, end='\r')


def log(*args, **kwargs):
    ts = datetime.now().strftime(ISO_FORMAT)[:-3]
    if 'ts' not in kwargs or kwargs['ts'] is not False:
        args = [ts, *args]
    if 'ts' in kwargs:
        del kwargs['ts']
    print(*args, **kwargs)
    if log_f:
        print(*args, **kwargs, file=log_f)
        log_f.flush()


def save_np_data(path, *args, **kwargs):
    if len(kwargs) == 0:
        kwargs = {'data': args[0]}
    np.savez_compressed(path + '.npz', **kwargs)


def load_np_data(path, key='data'):
    return np.load(path + '.npz')[key]


def save_csv(path, rows):
    with open(path, 'w') as f:
        csvWriter = csv.writer(f)
        csvWriter.writerows(rows)


def print_args(args):
    arg_dict = vars(args)
    m = max(len(k) for k in arg_dict.keys())
    log('Running args:')
    for k, v in arg_dict.items():
        print(f'  {k} {" " * (m - len(k))}: {v}')


def simplify_pat_ids(data):
    return ['-'.join(p.split('-')[:3]) for p in data]


def get_safe_path_obj(file_path):
    ftype = type(file_path)
    if not isinstance(file_path, Path):
        if ftype == str:
            return Path(file_path)
        else:
            raise TypeError(f'Type of file_path must be either str or pathlib.Path but given {ftype}')
    return file_path
