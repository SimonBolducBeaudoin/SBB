#!/bin/env/python
#! -*- coding: utf-8 -*-

import itertools
import numpy

def mklist(v):
    """
    Force v to be compatible with a list (either list, tuple or numpy.array)
    Similar to np.asarray but for list
    """
    v = v if isinstance(v, (list, tuple, numpy.ndarray)) else [v]
    return v
    
def super_enumerate(*args):
    """
        A multi-D equivalent of enumerate
        Iterates over the outer product of the args while
        
        Inputs : a list of 1D array
        
        Returns : 
            .next() returns a tuple of 2 tuples
            1st tuple  (...,k,j,i) of index    in the tensor representation 
            2nd tuple  (...,k,j,i) of elements in the tensor representation
    """
    if len(args) == 0 : # called empty
        return iter(()) , iter(())
    index_vec = ()
    for a in args :
        l = len(a)
        index_vec += ( range(l) , )
        
    it_idx = itertools.product(*index_vec) 
    it = itertools.product(*args)
    return zip( *(it_idx,it) )