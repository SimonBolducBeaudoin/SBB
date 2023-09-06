#!/bin/env/python
#! -*- coding: utf-8 -*-

"""Regroups function used to combine and clean up data produced by experiments"""

import numpy as _np

from SBB.Python_extra.python_extra import super_enumerate as _super_enumerate
from SBB.Numpy_extra.numpy_extra import argwhere_tuple as argwhere_tuple , find_where_A_equal_B as _find_AB

def get_untreated_files(folder,pattern= 'exp_??????-??????.npz',ReDo=False,AddIgnore=True,LookForLonelyIgnore=True): 
    """ 
    Used to get a list of files that have not been treated yet (aka not flagged with ignore yet). 
         
    Basic use (ReDo=False,AddIgnore=False,LookForLonelyIgnore=False) 
        Looks inside folder for files matching "pattern",  
        also looks for files matching "pattern"+.ignore  
        and return a list of pattern files excluding those that have a corresponding .ignore sister 
        Exemple :  
        If folder content is : 
             'exp_000000-000000.npz' 
             'exp_000000-000001.npz' 
             'exp_000000-000001.npz.ignore' 
        then it returns 
             ['exp_000000-000000.npz'] 
              
    Options : 
        ReDo (False by default) 
            Doesn't filter for .ignore file and hence returns all files matching the pattern 
        AddIgnore (True by default) 
            Adds .ignore to files name (modifies files names) 
        LookForLonelyIgnore 
            Looks if some "pattern.ignore" files dont have a matching "parttern" file and add them to the list 
    """ 
    pattern_list = glob.glob(folder+os.sep+pattern) 
    ignore_list  = glob.glob(folder+os.sep+pattern+'.ignore') 
     
    # New_list with only experiments that have not already been processed     
    not_done_list = [ exp for exp in pattern_list if (exp + '.ignore' not in ignore_list) ] 
 
    # Adds .ignore to the latter 
    if AddIgnore : 
        for exp in not_done_list : 
            os.rename(exp,exp+'.ignore') 
 
    # Get the "parttern.ignore" files that don't have a corresponding "pattern" file 
    lonely_ignore_list = [exp.replace('.ignore', '') for exp in ignore_list if ( (exp.replace('.ignore', '') not in pattern_list) and (LookForLonelyIgnore) ) ] 
     
    if ReDo : 
        return pattern_list + lonely_ignore_list 
    else : 
        return not_done_list + lonely_ignore_list 

def remove_zeros_subarrays(arr):
    """
    Returns only non-zeros elements along the firts dimension/axis = 0
    
    Example :
        Say arr [20,...] = 0.0 is full of zeros
        Then remove_zeros_subarrays(arr) will return arr with line 20 removed.
    """
    # Get the indices of all subarrays that do not contain zeros
    non_zero_indices = _np.argwhere(arr)
    non_zero_indices = _np.unique(non_zero_indices[:, 0])

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
    non_nan_indices = _np.argwhere(~_np.isnan(arr))
    non_nan_indices = _np.unique(non_nan_indices[:, 0])

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
        w = _argwhere_tuple(V,v)
        if w : # At least one match found
            indx = w[0] # Take only the first match's coordinate
            flat_indx = _np.ravel_multi_index(indx, V_shape) # convert to flat repr
            S_reps += [ S[:, flat_indx:flat_indx+1, ...] ,]
        else :
            S_reps += [ _np.nan ,]# place holder
    return S_reps
    
def sort_A_acording_to_B(A_list,B_list):
    """
    Sorts A_list according to the sorted values of B_list.
    Assumes A_list and B_list have identical shapes.

    Steps:
    1. Get the number of dimensions of the conditions in B_list.
    2. Sort the list of conditions in B_list in ascending order and keep track of the index of the first 
       occurrence of each unique value for each dimension.
    3. Get the sorted unique values of B for each dimension.
    4. Concatenate the values in A_list for each dimension and sort the resulting array based on how 
       B was sorted in step 2. Use the sorted indices from step 2 to sort A_list based on the sorted B.
    5. Return the sorted A and B.

    Arguments:
    A_list -- list of arrays to be sorted
    B_list -- list of arrays to be used as the sort keys

    Returns:
    A -- sorted A_list
    B -- sorted B_list
    """
    n_dim_cdn = len(B_list[0])
    # Sorts B (idx = 0) and keeps index of the first occurences (idx = 1)
    S = [ _np.unique(_np.concatenate([v[n] for v in B_list]),return_index=True) for n in range(n_dim_cdn) ]
    # get the sorted B for each dimension
    B   = [ s[0] for s in S ] 
    A   = [ (_np.concatenate([v[n] for v in A_list])[s[1]] ) for n,s in zip(range(n_dim_cdn),S) ]
    return A,B
    
def combine_arrays(S_list,V_list,remove_nans=True,remove_zeros=False):
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
    if isinstance( V_list[0][0], (list, tuple, _np.ndarray) ) :
        pass # List of list ok
    elif isinstance( V_list[0], (list, tuple, _np.ndarray) ) :
        # Single list of conditions
        V_list = [ [v] for v in V_list ]
    else :
        raise Exception("V_list has to be like a list of list of 1D np.arrays")
        
    # Check that Vs all have the same length
    if all( [ len(V_list[0]) != len(V) for V in V_list ] ) :
        raise Exception("len (V_0) != len(V_n)")
    
    if remove_nans :
        S_list = [ remove_nan_subarrays(Sn) for Sn in S_list ]
    if remove_zeros :
        S_list = [ remove_zeros_subarrays(Sn) for Sn in S_list ]
    
    # List of the number of repetitions for each Sn
    n_reps = [ len(Sn) for Sn in S_list ]
    
    # Get the number of dimensions of the experimental conditions
    n_dim_cdn = len(V_list[0])
    
    # A list of unique conditions for each experimental variables
    V =[ _np.unique(_np.concatenate([v[n] for v in V_list])) for n in range(n_dim_cdn) ]
    
    # For each experiment repetitions we can find the corresponding index in the final cdn array for each cdn axis
    # Index lists
    Idx_list = [[ _find_AB(v_cdn,v_unique)for v_cdn,v_unique in zip(V_rep,V)] for V_rep in V_list]
    
    # Shape of the unique combined experimental variables
    l = tuple( len(v) for v in V )
    
    # A list of list of lenght for each experimental variables of each experiments
    l_list = [[len(v) for v in V_rep] for V_rep in V_list ]
    
    # A List of reshaped measurment with all experimental variables flatten in a single dimension
    S_list = [ Sn.reshape( (Sn.shape[0],) + (_np.prod(ln),) + Sn.shape[1+n_dim_cdn:] ) for Sn,ln in zip(S_list,l_list) ]
    
    
    # The output array with the right dimension full of nans
    S_shape = (sum(n_reps), _np.prod(l), ) + S_list[0].shape[2:] # (n_reps, a,...,z            , ...                )
    S = _np.full(S_shape, _np.nan, dtype=S_list[0].dtype)        # (n_reps, ... conditions ... , ... other axis ... )
            
    n = 0
    for n_rep,s,idx_list in zip( n_reps,S_list,Idx_list) : # over inputs 
        for i, (from_idx, goto_idx) in enumerate( _super_enumerate(*idx_list) ):
            # Copy repetitions into the final combined array  
            j = _np.ravel_multi_index(goto_idx,l) 
            S[n:n+n_rep,j:j+1,...] = s[:,i:i+1,...]
        n += n_rep
                
        # Restoring shapes
    S.shape = (S.shape[0],) + l + S.shape[2:]
    
    # Return the combined arrays
    return V, S , n_reps
    
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
    Out = _np.full(C.shape,_np.nan,dtype=C.dtype)
    for i,t in enumerate(_np.isnan(C)):
        for j,tt in enumerate(t):
            if all(tt==False):
                break
        Out[i,:len(t)-j,...] = C[i,j:,...]
    return Out.swapaxes(0,axis)
    
#if __name__ == '__main__':
    # Scipt that test/show the functionnality of the library