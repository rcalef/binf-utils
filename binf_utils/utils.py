#! /usr/bin/env python3.3
import gzip
import bz2

def open_file(filename,mode="r"):
    if filename.endswith(".gz"): 
        return gzip.open(filename, mode)
    elif filename.endswith(".bz2"):
        return bz2.open(filename,mode)
    else:
        return open(filename,mode)
