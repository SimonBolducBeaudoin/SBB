#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np

from SBB.Python_extra.python_extra import mklist as _mklist, super_enumerate as _super_enumerate

def get_index(Xs,x):
    """
    Returns the single index for wich Vdc = V
    Returns false if it doesn't exist
    """
    tmp = where(Vdc==V)[0]
    return False if tmp.size==0 else tmp[0]

def find_nearest_A_to_a(a,A):
    a = _np.array([a]) if type(a)==float else _np.array(a)
    A = _np.array(A)
    a_shape = a.shape
    a       = a.flatten()
    X       = _np.empty(a.shape)
    X_idxs  = _np.empty(a.shape,dtype=int)
    for i,x in enumerate(a) :
        index       = _np.abs(A-x).argmin()
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
    return _np.concatenate((X[...,-1:0:-1],X),axis=-1)

def cyclic_tansformation(X):
    """
        Cyclic translation of nd array
    """
    shape     = X.shape
    out       = _np.zeros(X.size)
    flat_input = X.flatten()
    out[0:-1] = flat_input[1:]
    out[-1]   = flat_input[0]
    out.shape = shape
    return out

def build_array_of_objects(shape,constructor,*args,**kargs):
    A = _np.r_[[constructor(*args,**kargs) for i in range(_np.prod(shape))]]
    A.shape = shape
    return A

class sub_flatiter(object) :
    """
    Iterates over a flat representation of arr up to axis.
    
    Ex :
    Say arr.shape = (2,3,4,5)
    
    it = sub_flatiter(arr,axis=-2)
    iterates over the first dimension of an array of shape (6,4,5)
    without modifying the initial array.
    """
    def __init__(self,arr,axis=-1):
        self.arr        = arr
        self.axis       = axis
        self.flat_shape = (int(_np.prod(arr.shape[:axis])),) + arr.shape[axis:]
        self.j          = -1
        self.j_max      = self.flat_shape[0] - 1
    def __iter__(self):
        return self
    def __next__(self):
        if (self.j >= self.j_max -1) :
            raise StopIteration
        self.j +=1
        # Flattening the arrays without copying.
        shape     = sub_flatten_no_copy(self.arr,axis=self.axis) 
        out       = self.arr[self.j,...]
        # Restoring the array's shape
        self.arr.shape = shape      
        return out 
        
def sub_flat(arr,axis=-1):
    """
    Returns an iterator that iterates over a flat representation of arr up to axis.
    
    Similar to flat
    
    Ex :
    Say arr.shape = (2,3,4,5)
    
    for a in sub_flat(arr,axis=-2):
        iterates over the first dimension of an array of shape (6,4,5)
        without modifying the initial array.
    """
    return sub_flatiter(arr,)

def sub_flatten(arr,axis=-1):
    return arr.reshape( (int(_np.prod(arr.shape[:axis])),) + arr.shape[axis:] ) 

def sub_flatten_no_copy(arr,axis=-1):
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
    arr.shape = (int(_np.prod(shape[:axis])),) + shape[axis:]
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
        max_shape = _np.maximum(max_shape, array(S.shape))

    # Create a padding tuple with the same number of dimensions as S
    if axis == None : # Pad on all axis
        padding = [ [(0, max_shape[i] - S.shape[i]) for i in range(len(S.shape))] for S in S_reps ]
    else : # Pad on specified axis
        axis = _mklist(axis)
        padding = [ [ (0, max_shape[i] - S.shape[i]) if i in axis else (0,0) for i in range(len(S.shape)) ] for S in S_reps ]
    
    # Pad all the input arrays with Nans so that they have the same shape
    padded_arrays = [ _np.pad(S, p, mode="constant", constant_values=float("nan")) for S,p in zip(S_reps,padding) ]

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
    return [v_idx for v_idx , v0  in _super_enumerate(*V) if v0 == v]
    
def where_tuple(V,v):
    """
    Same as argwhere but the output is grouped by dimensions and not by element, similarly to the output of np.where
    """
    # list of position tuples where the conditions is matched
    w = [v_idx for v_idx , v0  in _super_enumerate(*V) if v0 == v] 
    # converting to right output format    
    return [ [x[dim]    for x in w] for dim in range(len(V))]
    
def find_where_A_equal_B(A,B):
    """
    Return the index of B where B[i] == A for each element of A.
    Intended for 1D arrays only
    """
    return [i for a in  A for i,b in enumerate(B) if a == b]
     
#if __name__ == '__main__':
    # Scipt that test/show the functionnality of the library