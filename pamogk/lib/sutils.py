import argparse
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


def ensure_suffix(path, suffix='.npz'):
    if isinstance(path, Path):
        return path.with_suffix(suffix)
    if isinstance(path, str):
        if path[-len(suffix):] != suffix:
            return path + suffix
        return path
    raise TypeError('path must be either pathlib.Path or str')


def np_save_npz(path, data=None, **kwargs):
    if len(kwargs) == 0:
        if data is None:
            raise ValueError('Either data or kwargs must be given')
        kwargs = {'data': data}
    path = ensure_suffix(path)
    np.savez_compressed(path, **kwargs)


def np_load_data(path, key=None):
    """
    Loads numpy data from npz files
    Parameters
    ----------
    path
    key: str 'data'
    If given loads data in that key, if None is given loads all of the data in npz

    Returns
    -------

    """
    data = np.load(path)
    if key is not None:
        if key not in data:
            log(f'Given key={key} is not in loaded data on path={path}')
        return data[key]
    return data


def save_csv(path, rows):
    with open(path, 'w') as f:
        csvWriter = csv.writer(f)
        csvWriter.writerows(rows)


def print_args(args, ignore_keys=None):
    if ignore_keys is None:
        ignore_keys = ['self']
    arg_dict = args
    if type(args) is argparse.Namespace:
        arg_dict = vars(args)
    arg_dict = dict((k, v) for k, v in arg_dict.items() if k not in ignore_keys)
    m = max(len(k) for k in arg_dict.keys())
    log('Running args:')
    for k, v in arg_dict.items():
        log(f'  {k} {" " * (m - len(k))}: {v}')


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


def str2bool(v):
    """
    This is useful for the case of giving boolean parameters
    Parameters
    ----------
    v

    Returns
    -------

    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def str2path(v):
    """
    Just a wrapper for
    Parameters
    ----------
    v

    Returns
    -------

    """
    try:
        return Path(v)
    except TypeError:
        raise argparse.ArgumentTypeError('Invalid path provided')
