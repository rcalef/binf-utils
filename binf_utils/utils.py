import gzip
import bz2
import sys

def open_file(filename,mode="r"):
    if filename.endswith(".gz"): 
        return gzip.open(filename, mode)
    elif filename.endswith(".bz2"):
        return bz2.open(filename,mode)
    elif '-' == filename:
        if mode == 'r':
            return sys.stdin
        elif mode == 'w':
            return sys.stdout
    else:
        return open(filename,mode)

def max_argmax(**kwargs):
    arg = max(kwargs.keys(), key = lambda x: kwargs[x])
    return (kwargs[arg],arg)
