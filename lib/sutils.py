import time
import os
import sys

class timeit(object):
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

    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        print("Started:", self.f.__name__)
        t = time.time()
        res = self.f(*args, **kwargs)
        print("Finished: {} elapsed: {:.2f}s".format(self.f.__name__, time.time() - t))
        return res

def safe_create_dir(d):
    if not os.path.exists(d):
        print('Dir not found creating:', d)
        os.makedirs(d)

# _print = print

# def print(*args, flush=True, **kwargs):
#     _print(*args, **kwargs)
#     if flush: sys.stdout.flush()
