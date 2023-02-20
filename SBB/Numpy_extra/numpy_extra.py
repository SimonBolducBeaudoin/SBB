#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy

from SBB.Python_extra.python_extra import mklist,super_enumerate

def get_index(Xs,x):
    """
    Returns the single index for wich Vdc = V
    Returns false if it doesn't exist
    """
    tmp = where(Vdc==V)[0]
    return False if tmp.size==0 else tmp[0]

def find_nearest_A_to_a(a,A):
    a = numpy.array([a]) if type(a)==float else numpy.array(a)
    A = numpy.array(A)
    a_shape = a.shape
    a       = a.flatten()
    X       = numpy.empty(a.shape)
    X_idxs  = numpy.empty(a.shape,dtype=int)
    for i,x in enumerate(a) :
        index       = numpy.abs(A-x).argmin()
        X_idxs[i]   = index 
        X[i]        = A[index]
    X.shape         = a_shape
    X_idxs.shape    = a_shape
    return X , X_idxs 

def symetrize(X):
    """
    Symetrize data along the last axis assuming the first point must not be duplicated.
    
    Input
    -----
        S2 : np.array with shape = (...,N)
    Output
        np.array with shape = (...,2*N-1)  
    """
    return numpy.concatenate((X[...,-1:0:-1],X),axis=-1)

def cyclic_tansformation(X):
    """
        Cyclic translation of nd array
    """
    shape     = X.shape
    out       = numpy.zeros(X.size)
    flat_input = X.flatten()
    out[0:-1] = flat_input[1:]
    out[-1]   = flat_input[0]
    out.shape = shape
    return out

def build_array_of_objects(shape,constructor,*args,**kargs):
    A = numpy.r_[[constructor(*args,**kargs) for i in range(numpy.prod(shape))]]
    A.shape = shape
    return A

def sub_flatten(arr,axis=-1):
    """
    Equivalent to 
    arr = arr.reshape( (int(numpy.prod(arr.shape[:axis])),) + arr.shape[axis:] ) 
    but doesn't make a copy of the array instead it modifies arr.shape.
    
    Warning ! 
    ---------
    This function will modify the arr input's shape. 
    
    Returns
    ------- 
    old_shape : tuple
        Can be used to restore the initial shape of the array afterward.
        It it the responsability of the user to do so.
    """
    
    shape    = arr.shape
    arr.shape = (int(numpy.prod(shape[:axis])),) + shape[axis:]
    return shape 
    
def pad_arrays( *S_reps, **kwargs ):
    """
    Pads all arrays with Nans such that they all have the shape of the bigest array.
    
    Parameters
    ----------
    S_reps : list of mp.array
    
    axis : (falcultative) list of int
        The axis along which to do enlarge the arrays
    
    Returns
    --------
    padded_arrays : list of np.array of the same shape
    """
    axis = kwargs.get('axis')
    # Find the maximum shape of all the input arrays
    max_shape = array(S_reps[0].shape)
    for S in S_reps[1:]:
        max_shape = numpy.maximum(max_shape, array(S.shape))

    # Create a padding tuple with the same number of dimensions as S
    if axis == None : # Pad on all axis
        padding = [ [(0, max_shape[i] - S.shape[i]) for i in range(len(S.shape))] for S in S_reps ]
    else : # Pad on specified axis
        axis = mklist(axis)
        padding = [ [ (0, max_shape[i] - S.shape[i]) if i in axis else (0,0) for i in range(len(S.shape)) ] for S in S_reps ]
    
    # Pad all the input arrays with Nans so that they have the same shape
    padded_arrays = [ numpy.pad(S, p, mode="constant", constant_values=float("nan")) for S,p in zip(S_reps,padding) ]

    return padded_arrays
    
def argwhere_tuple(V,v):
    """
    Find where V == v in the outer product tuple space.
    The strength of this approach is that works event if arrays V[i] and V[j] have different type not compatible for comparisons together.
    
    Similarly to argwhere matches are grouped by element (see example).
    
    Ex :
        >>> V = [['A', 'B'], array([20, 22]), [float, float]]
        >>> v = ('A', 22, float)
        >>> indxs = argwhere_tuple(V,v)
        >>> indxs[0] # The coordinates in the outter product space for the first match
        (0, 1, 0)
        >>> indxs[1] # The coordinates in the outter product space for the 2nd match
        (0, 1, 1)
    Parameters
    ----------
    V : a list of 1D array-like
    v : a tuple with len(v) = len(V)
    
    Returns
    --------
    inxs : list 
        A list of match where outproduct of V == v
    """
    return [v_idx for v_idx , v0  in super_enumerate(*V) if v0 == v]
    
def where_tuple(V,v):
    """
    Same as argwhere but the output is grouped by dimensions and not by element, similarly to the output of np.where
    """
    # list of position tuples where the conditions is matched
    w = [v_idx for v_idx , v0  in super_enumerate(*V) if v0 == v] 
    # converting to right output format    
    return [ [x[dim]    for x in w] for dim in range(len(V))]
     
     

       

#if __name__ == '__main__':
    # Scipt that test/show the functionnality of the library