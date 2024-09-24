#!/bin/env/python
#! -*- coding: utf-8 -*-

"""
This saves results in a dictionnary in ram and re-use the result if the same call is made again.

I could add a similar behevior for heavy function but saving to disk in hte form of a np array perhaps ?
"""

import sys as _sys
__all__ = ["MemoizeMutable"]

import pickle
class MemoizeMutable:
    def __init__(self, fn, verbose=False):
        self.fn = fn
        self.memo = {}
        self.verbose=verbose
    def __call__(self, *args, **kwds):
        str = pickle.dumps(args, 1)+pickle.dumps(kwds, 1)
        if not str in self.memo: 
            self.memo[str] = self.fn(*args, **kwds)
        return self.memo[str] 

