#!/bin/python
# -*- coding: utf-8 -*-

import sys, os, platform, time
import numpy as np
import matplotlib.pyplot as plt
from numpy import uint8, int8, uint16, int16, double
from numpy import ndarray, ceil, log2, iinfo, zeros, allclose, arange, array 
from numpy import floor, log10, savez_compressed, load
from decimal import Decimal


# Uncomment if you're not installing it in SBB library
#import acorrs_wrapper
#from acorrs_wrapper import set_mpreal_precision

# For compatibility with my own installation
import SBB.AutoCorr.acorrs_wrapper as acorrs_wrapper
from SBB.AutoCorr.acorrs_wrapper import set_mpreal_precision

# Applies to instances created afterwards
set_mpreal_precision(48)

# For automatic fftchunk determination
def closest_power_of_two(x):
    a = 2**int(ceil(log2(x)))
    b = 2**int(floor(log2(x)))
    if abs(a-x)<abs(x-b):   # Biased towards smallest value if dead center
        return a
    else:
        return b

# Returns the proper class. Fancy name: factory. Ghetto name: wrapper wrapper.
def ACorrUpTo(k, data, phi=False, fft=None, fftchunk='auto', k_fft=32, k_fft_factor=16):
    if type(data) is ndarray:
        dtype = data.dtype.name
    else:
        dtype = data
    
    if phi:
        fft = False

    if fft is None: 
        if k>=k_fft:  # k_fft is empirical for each system 
            fft = True
        else:
            fft = False

    if fftchunk == 'auto':
        # Since fftwlen is a power of 2, fftchunk is optimal if it is too
        fftchunk = closest_power_of_two(k_fft_factor*k)

    if fft and k>fftchunk:
        fftchunk = int(2**ceil(log2(k))) # Ceil to power of two
    
    classname = "ACorrUpTo{fft}_{dtype}".format(dtype=dtype, fft="FFT" if fft else "Phi" if phi else "")
    
    if fft:
        retClass = getattr(acorrs_wrapper, classname)(k, fftchunk)
    elif phi:
        retClass = getattr(acorrs_wrapper, classname)(k, phi)
    else:
        retClass = getattr(acorrs_wrapper, classname)(k)
    
    if type(data) is ndarray:
        retClass(data)
    
    return retClass


# For testing

# Computes phase-resolved a.res using Decimal accumulators
# Casting result to double should be exactly a.res[k]
def check_ak(a,k):
    nk = a.n-k
    rk = a.rk[k]
    m = a.m
    bk = a.bk[k]
    gk = a.gk[k]
    return (rk-(m-bk)*(m-gk)/nk)/nk

# Computes phase-resolved a.res using Decimal accumulators
# Casting result to double should be exactly a.res[f,k]
def check_afk_phi(a,f,k):
    nfk = a.nfk[f,k]
    rfk = a.rfk[f,k]
    mf = a.mf[f]
    mfpk = a.mf[(f+k)%a.l]
    bfk = a.bfk[f,k]
    gfk = a.gfk[f,k]
    return (rfk-((mf-bfk)*(mfpk-gfk))/nfk)/nfk

# Computes phase-resolved a.res0 using Decimal accumulators
# Casting result to double should be exactly a.res0[k]
def check_ak_phi(a,k):
    nfk = a.nfk.sum(axis=0)[k]
    rfk = a.rfk.sum(axis=0)[k]
    mf = a.mf.sum()
    mfpk = a.mf.sum()
    bfk = a.bfk.sum(axis=0)[k]
    gfk = a.gfk.sum(axis=0)[k]
    return (rfk-((mf-bfk)*(mfpk-gfk))/nfk)/nfk

# Converts an ACorrUpTo object to a dict with the same information
def a_to_dict_phi(a):
    ks = 'bk block_processed chunk_processed chunk_size gk k m n res rk'.split(' ')
    return {k:getattr(a,k) for k in ks}
    
# Converts an ACorrUpToFFT object to a dict with the same information
def a_to_dict_fft(a):
    ks = 'bk block_processed chunk_processed chunk_size counter_max fftwlen gk k len m n res rk'.split(' ')
    return {k:getattr(a,k) for k in ks}    

# Converts an ACorrUpToPhi object to a dict with the same information
def a_to_dict_phi(a):
    ks = 'bfk bk block_processed chunk_processed chunk_size gfk gk k l mf n nfk res res0 rfk'.split(' ')
    return {k:getattr(a,k) for k in ks}

    
