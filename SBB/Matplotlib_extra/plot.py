#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np

def nanplot(ax,x,y,**kwargs):
    not_nan = ~(_np.isnan(y)) 
    return ax.plot(x[not_nan],y[not_nan],**kwargs)

def plot_interval(ax,x,Y,std=True,n_sigma=3,min_max=True,label=None,linestyle=None,marker=None,alpha=[None,0.75,0.5], **kwargs):
    """
    Almost the same as plot but plota Y.mean(axis=0) and (optionnaly) some visual metrics on deviation.
    
    Inputs
    ------
    ax : an axis
    x  : the abscisse
    Y  : np.array of shape (N,len(x)) 
    
    std : plot the std
    n_sigma : the number of std deviation to be used
    min_max : plot the min_max
    
    **kwargs will be passed to all subplots
    label and linestyles applies to the average plot
    alpha controls the transparency for deviation plots
    
    returns
    -------
    line, same output as ax.plot
    """
    Y_m = _np.nanmean( Y,axis=0 )
    not_nan = ~(_np.isnan(Y_m)) 
    
    all_nan = all( not_nan ==False)
    
    if all_nan :
        line = None
    else :
        line, = ax.plot(x[not_nan], Y_m[not_nan] ,label=label,marker=marker,linestyle=linestyle,**kwargs)
        if std :
            Y_std = _np.nanstd( Y,axis=0 )
            ax.fill_between(x[not_nan], (Y_m + n_sigma*Y_std)[not_nan], (Y_m - n_sigma*Y_std)[not_nan], facecolor = line.get_color(),alpha=alpha[1],**kwargs)
        if min_max :    
            Y_max = _np.nanmax( Y,axis=0 )
            Y_min = _np.nanmin( Y,axis=0 )
            ax.fill_between(x[not_nan], Y_max[not_nan],Y_min[not_nan], facecolor = line.get_color(),alpha=alpha[2],**kwargs)
    return line, 