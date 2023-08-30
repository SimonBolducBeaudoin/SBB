#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy
from SBB.Numpy_extra.numpy_extra import sub_flatten_no_copy

"""
This module is intended to help with the use of SBB.Histogram.histograms

Last update
-----------
    added compute_moments
Todos : 
Bugs :
"""

def compute_moments_par(Hs,x,order = 8,Cxs=None):
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
        
    def _moments(h,i,j,bins,order):
        moments[i,j,:] = h.std_moments(bins,order,no_clip=True)
        
    moments = numpy.full( Hs.shape + (order+1,), numpy.nan ) 
    ms_shape = sub_flatten_no_copy(moments,axis=-2) # this function modifies the input's shape
    Hs_shape = sub_flatten_no_copy(Hs,axis=-1)      # this function modifies the input's shape
    
    threads = [ threading.Thread(target= _moments , args = (h,i,j,x*Cx,order)) for i, (H,Cx) in enumerate( zip( Hs, Cxs.flat ) ) for j,h in enumerate(H)]        
    for th in threads :
        th.start()
    for th in threads :
        th.join()
    #Restore shapes
    Hs.shape = Hs_shape
    moments.shape = ms_shape
    return moments
    
# def compute_moments_par(Hs,x,order = 8,Cxs=None):
    # """
    # Computes standardized moments directly from Histogram class
    # see: SBB.Histogram.moments_cumulants import std_moments 
    
    # Inputs
    # ------
    # Hs : numpy array of 1D histograms objects
        # Hs.shape  = (..., n )
    # x  : numpy array of the center for each bin (see: SBB.Histogram.Histogram_uint64_t_double.abscisse)
        # x.shape = (n,)
    # Cxs : numpy array of corrections the abscisse for each kernel
        # used to help with half normalisation
        # Cxs.shape = Hs.shape[:-1]
    # order : computing standardized moments up to this order
    # """ 
    # if not(Cxs):
        # Cxs = numpy.ones( Hs.shape[:-1] )
        
    # def _moments(h,i,j,bins,order):
        # moments[i,j,:] = h.std_moments(bins,order,no_clip=True)
        
    # moments = numpy.full( Hs.shape + (order+1,), numpy.nan ) 
    # ms_shape = sub_flatten(moments,axis=-2) # this function modifies the input's shape
    # Hs_shape = sub_flatten(Hs,axis=-1)      # this function modifies the input's shape
    # for i, (H,Cx) in enumerate( zip( Hs, Cxs.flat ) ): 
        # bins = x*Cx
        # threads = [ threading.Thread(target= _moments , args = (h,i,j,bins,order)) for j,h in enumerate(H)]
        # for th in threads :
            # th.start()
        # for th in threads :
            # th.join()
    ##Restore shapes
    # Hs.shape = Hs_shape
    # moments.shape = ms_shape
    # return moments
    
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
    if Cxs is None :
        Cxs = numpy.ones( Hs.shape[:-1] )
        
    moments = numpy.full( Hs.shape + (order+1,), numpy.nan ) 
    ms_shape = sub_flatten_no_copy(moments,axis=-2) # this function modifies the input's shape
    Hs_shape = sub_flatten_no_copy(Hs,axis=-1)      # this function modifies the input's shape
    for i, (H,Cx) in enumerate( zip( Hs, Cxs.flat ) ): 
        bins = x*Cx
        for j, h  in enumerate( H ) : 
            # Computes standardized moments up to order 
            moments[i,j,:] = h.std_moments(bins,order,no_clip=True)
    #Restore shapes
    Hs.shape = Hs_shape
    moments.shape = ms_shape
    return moments
   
def std_moments_to_moments(moments):
    for i in range(moments.shape[-1]):
        if i > 2 :
            moments[...,i] *= moments[...,2]**(float(i)/2.0)
    return moments