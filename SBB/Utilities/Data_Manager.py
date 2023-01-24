#!/bin/env/python
#! -*- coding: utf-8 -*-

"""Regroups function used to combine and clean up data produced by experiments"""

import numpy

from SBB.Utilities.python_extra import super_enumerate
from SBB.Utilities.numpy_extra import argwhere_tuple

def remove_zeros_subarrays(arr):
    """
    Returns only non-zeros elements along the firts dimension/axis = 0
    
    Example :
        Say arr [20,...] = 0.0 is full of zeros
        Then remove_zeros_subarrays(arr) will return arr with line 20 removed.
    """
    # Get the indices of all subarrays that do not contain zeros
    non_zero_indices = numpy.argwhere(arr)
    non_zero_indices = numpy.unique(non_zero_indices[:, 0])

    # Return the subarrays at the non-Nan indices
    return arr[non_zero_indices]
    
def remove_nan_subarrays(arr):
    """
    Returns only non-nan elements along the firts dimension/axis = 0
    
    Example :
        Say arr [20,...] = nan is full of nan
        Then remove_nan_subarrays(arr) will return arr with line 20 removed.
    """
    # Get the indices of all subarrays that do not contain Nans
    non_nan_indices = numpy.argwhere(~isnan(arr))
    non_nan_indices = numpy.unique(non_nan_indices[:, 0])

    # Return the subarrays at the non-Nan indices
    return arr[non_nan_indices]
    
def repetitions_of_v(S_list,V_list,v):
    """
    Returns a list, S_reps, of experimetnal repetitions of the condition v
    
    S_reps[i] have shape (:,1,...) or are simply nan is no match was found.
    
    Notes
    -----
    We assume that a given v apears only once per experiment
    """
    S_reps = []
    # A liste of where( vi == v) for each Vi
    for S,V in zip(S_list,V_list) :
        V_shape = [len(vv) for vv in V]
        w = argwhere_tuple(V,v)
        if w : # At least one match found
            indx = w[0] # Take only the first match's coordinate
            flat_indx = numpy.ravel_multi_index(indx, V_shape) # convert to flat repr
            S_reps += [ S[:, flat_indx:flat_indx+1, ...] ,]
        else :
            S_reps += [ numpy.nan ,]# place holder
    return S_reps
    
def combine_arrays(S_list,V_list):
    """
    Combines experimental data (S_list) in order not to loose repetitions of experimental conditions.
    
    Parameters
    ----------
    S_list is a list of data each element coming from a different experiment.
        S_list[k].shape = (n_reps,n_cdn_0,...,n_cdn_i,...)
    V_list is a list of list containing each experimental axis
        len( V_list[k] ) = i the number of experimental axis
        len( V_list[0][0] ) = n_cdn_0_0
        ...
        len( V_list[0][1] ) = n_cdn_0_1
        ...
        len( V_list[1][0] ) = n_cdn_1_0
        ...
            
    Returns
    -------
    V: list
        List of the combined experimental conditions.
    S: numpy array
        Shape (n_rep, combined_cdn_0, ..., combined_cdn_i , ...)
        
    Notes
    ------
    A giving axis of coditions V_list[i] must only containt the same condition once.
    The behaviour is not guarantteed if this is not respected.    
    """     
    # Check that V_list like a list of list of 1D arrays
    if isinstance( V_list[0][0], (list, tuple, numpy.ndarray) ) :
        pass # List of list ok
    elif isinstance( V_list[0], (list, tuple, numpy.ndarray) ) :
        # Single list of conditions
        V_list = [ [v] for v in V_list ]
    else :
        raise Exception("V_list has to be like a list of list of 1D np.arrays")
        
    # Check that Vs all have the same length
    if all( [ len(V_list[0]) != len(V) for V in V_list ] ) :
        raise Exception("len (V_0) != len(V_n)")
    
    # Get the number of dimensions of the experimental conditions
    n_dim_cdn = len(V_list[0])
    
    # A list of unique conditions for each experimental variables
    V =[ numpy.unique(numpy.concatenate([v[n] for v in V_list])) for n in range(n_dim_cdn) ]
    # Shape of the unique combined experimental variables
    l = tuple( len(v) for v in V )
    
    # A list of list of lenght for each experimental variables of each experiments
    l_list = [[len(v) for v in V_rep] for V_rep in V_list ]
    
    # A List of reshaped measurment with all experimental variables flatten in a single dimension
    S_list = [ Sn.reshape( (Sn.shape[0],) + (numpy.prod(ln),) + Sn.shape[1+n_dim_cdn:] ) for Sn,ln in zip(S_list,l_list) ]

    # List of the number of repetitions
    n_reps = [ len(Sn) for Sn in S_list ]
    
    # The output array with the right dimension full of nans
    S_shape = (sum(n_reps), numpy.prod(l), ) + S_list[0].shape[2:]
    S = numpy.full(S_shape, numpy.nan, dtype=S_list[0].dtype)
    
    # Iterate over the combined experimental conditions
    for i, (v_idx, v) in enumerate( super_enumerate(*V) ):
        # A list og repetitions of v 
        S_reps = repetitions_of_v(S_list,V_list,v)
        # Copy repetitions into the final combined array
        n = 0
        for n_rep, s_rep in zip ( n_reps , S_reps ) :
            S[n:n+n_rep,i:i+1,...] = s_rep
            n += n_rep
            
    # Restoring shapes
    S.shape = (S.shape[0],) + l + S.shape[2:]

    # Return the combined arrays
    return V, S
    
def nan_collapse(S,axis=0):
    """
    Removes leading nans along axis 
    
    axis = 0
    
    | nan 0   1  2  3         | 0 1 2  3   nan nan
    | nan nan 5  6  7     ==> | 5 6 7  nan nan nan
    | 8   9   10 11 23        | 8 9 10 11  23 
    
    """
    C   = S.copy()
    C   = C.swapaxes(0,axis)
    Out = numpy.full(C.shape,numpy.nan,dtype=C.dtype)
    for i,t in enumerate(numpy.isnan(C)):
        for j,tt in enumerate(t):
            if all(tt==False):
                break
        Out[i,:len(t)-j,...] = C[i,j:,...]
    return Out.swapaxes(0,axis)
    
#if __name__ == '__main__':
    # Scipt that test/show the functionnality of the library