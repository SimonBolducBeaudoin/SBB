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
    
    all_nan = _np.all( not_nan ==False)
    
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

def scope_interval_update(ax,x,Y,std=True, n_sigma=3, min_max=True, label=None, linestyle=None, marker=None, alpha=[None, 0.75, 0.5],line_idx=0) :
    col = ax.collections
    Y_m = _np.nanmean( Y,axis=0 )
    not_nan = ~(_np.isnan(Y_m))
    ax.lines[line_idx].set_data(x,Y_m)
    if std :
        Y_std = _np.nanstd( Y,axis=0 )
        verts = verts_to_fill_between(x[not_nan], (Y_m + n_sigma*Y_std)[not_nan], (Y_m - n_sigma*Y_std)[not_nan])
        col[line_idx*2  ].set_verts(verts[None,...])
    if min_max :    
        Y_max = _np.nanmax( Y,axis=0 )
        Y_min = _np.nanmin( Y,axis=0 )
        verts = verts_to_fill_between(x[not_nan], Y_max[not_nan], Y_min[not_nan]) 
        col[line_idx*2+1].set_verts(verts[None,...])

def verts_to_fill_between(x, y1, y2=0,closed=True):
    """
    Computes the vertices for fill_between 
    """
    x, y1, y2 = np.broadcast_arrays(np.atleast_1d(x), y1, y2)
    
    N = len(x)
    X = np.zeros((2 * N + 2, 2), float)
    
    start = x[0], y2[0]
    end = x[-1], y2[-1]

    X[0] = start
    X[N + 1] = end

    X[1:N + 1, 0] = x
    X[1:N + 1, 1] = y1
    X[N + 2:, 0] = x[::-1]
    X[N + 2:, 1] = y2[::-1]
    
    if closed : # closing the polygone
        X = np.concatenate([X, X[0:1]]) 
    return X  

def scope_update(ax,x,y,line_idx=0,option=None) :
    """
    Example of use :
    
    def plot_YvsX_then_refresh_auto(ax,x,Y,option='Stack'):
        if ax is None :
            fig,ax = subplots(1,1)
        if len(ax.lines)==0 :
            ax.plot(x,y)
            ax.set_ylim(...)
            ax.grid(True)
            ax.set_xlabel('x[~]')
            ax.set_ylabel('y[~]')
        else :    
            scope_update(x,y,ax,option)
        return ax
        
    """
    if option == 'Stack' : # adds a lines each time the function is called
        pc  =  ax._get_lines.prop_cycler
        line = matplotlib.lines.Line2D( x,y, **pc.next() )
        ax.add_line(line)
    elif option == 'keep_first' :
        if len(ax.lines) <= 1 :
            pc  =  ax._get_lines.prop_cycler
            line = matplotlib.lines.Line2D( x,y, **pc.next() )
            ax.add_line(line)
        else :
            ax.lines[1].set_data(x,y)
    else : # overwirte data
        ax.lines[line_idx].set_data(x,y)    