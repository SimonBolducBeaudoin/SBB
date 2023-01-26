#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy
from SBB.Utilities.numpy_extra import sub_flatten 

"""
This module is intended to help with the use of SBB.Histogram.histograms

Last update
-----------
    added compute_moments
Todos : 
Bugs :
"""

def compute_moments(Hs,x,order = 8,Cxs=None):
    """
    Computes standardized moments directly from Histogram class
    see: SBB.Histogram.moments_cumulants import std_moments 
    
    Inputs
    ------
    Hs : numpy array of 1D histograms objects
        Hs.shape  = (..., n )
    x  : numpy array of the center for each bin (see: SBB.Histogram.Histogram_uint64_t_double.abscisse)
        x.shape = (n,)
    Cxs : numpy array of corrections the abscisse for each kernel
        used to help with half normalisation
        Cxs.shape = Hs.shape[:-1]
    order : computing standardized moments up to this order
    """ 
    if not(Cxs):
        Cxs = numpy.ones( Hs.shape[:-1] )
        
    moments = numpy.full( Hs.shape + (order+1,), numpy.nan ) 
    ms_shape = sub_flatten(moments,axis=-2) # this function modifies the input's shape
    Hs_shape = sub_flatten(Hs,axis=-1)      # this function modifies the input's shape
    for i, (H,Cx) in enumerate( zip( Hs, Cxs.flat ) ): 
        bins = x*Cx
        for j, h  in enumerate( H ) : 
            # Computes standardized moments up to order 
            moments[i,j,:] = h.std_moments(bins,order,no_clip=True)
    #Restore shapes
    Hs.shape = Hs_shape
    moments.shape = ms_shape
    return moments