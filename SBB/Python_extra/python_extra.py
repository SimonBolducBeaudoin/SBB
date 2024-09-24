#!/bin/env/python
#! -*- coding: utf-8 -*-

import itertools as _ittools
import numpy as _np

def formated_tuple(frmt='{:0.2f}',T=tuple()):
    s = ''
    T = _np.array([T]) if (isinstance(T, float) or isinstance(T, _np.float64)) else T 
    for t in T :
        s += frmt.format(t)
        s += ', '
    return s
   
def mklist(v):
    """
    Force v to be compatible with a list (either list, tuple or numpy.array)
    Similar to np.asarray but for list
    """
    v = v if isinstance(v, (list, tuple, _np.ndarray)) else [v]
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
        
    it_idx = _ittools.product(*index_vec) 
    it = _ittools.product(*args)
    return zip( *(it_idx,it) )
    
def clean_function_source(func):
    """
    Returns the cleaned up source code of a function.

    Parameters:
    func (function): The function to get the source code for.

    Returns:
    str: The cleaned up source code of the function.
    """
    src = inspect.getsource(func)
    src = src[src.index(':\n') + 1:] # Removes everything before the first ":\n"
    src = re.sub(r'""".*?"""|\'\'\'.*?\'\'\'', '', src, flags=re.DOTALL) # Removes comments 
    src = re.sub(r'#.*', '', src) # Removes comments 
    src = re.sub(r'\n\s*\n', '\n', src) # Removes empty lines
    src = re.sub(r' (?<!\t) +(?!\t)', ' ', src) # Removes non-essential spaces
    return src

def compare_function_sources(src1, src2):
    """
    Compare the cleaned up source code of two functions.

    Parameters:
    src1 (str): The cleaned up source code of the first function.
    src2 (str): The cleaned up source code of the second function.
    Returns:
    list: A list of line numbers of any differences between the two source code strings. If the source code strings are the same, the list will be empty.
    """
    # Split the source code strings into lines
    lines1 = src1.split('\n')
    lines2 = src2.split('\n')

    # Initialize a list to store the line numbers of any differences
    diff_lines = []

    # Compare the lines of the source code
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        if line1 != line2:
            diff_lines.append(i + 1)  # Add the line number to the list

    # Return the list of line numbers
    return diff_lines    