#!/bin/env/python
#! -*- coding: utf-8 -*-

def get_twin(ax, axis='x'):
    """
    get_twin(ax, "x") returns ax's twinx if it has one.
    get_twin(ax, "y") returns ax's twiny if it has one.
    Source: 
    https:// stackoverflow.com /questions/ 36209575/ how-to-detect-if-a-twin-axis-has-been-generated-for-a-matplotlib-axis
    """
    assert axis in ("x", "y")        
    siblings = getattr(ax, "get_shared_"+axis+"_axes")().get_siblings(ax)
    for sibling in siblings:
        if sibling.bbox.bounds == ax.bbox.bounds and sibling is not ax:
            return sibling 
    return None
    
def get_twins(ax, axis='x'):
    """
    same as get_twin but return the whole familly of twins
    """
    assert axis in ("x", "y")        
    siblings = getattr(ax, "get_shared_"+axis+"_axes")().get_siblings(ax)
    twins = []
    for sibling in siblings:
        if sibling.bbox.bounds == ax.bbox.bounds and sibling is not ax:
            twins += [sibling,]
    return None if len(twins) == 0 else twins 
    
    
def get_next_ax(ax):
    """
    Returns the next ax in the current figure if it exist else returns None 
    """ 
    found_current_ax = False
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            found_current_ax = True
            continue
        if found_current_ax :
            return other_ax
    return None