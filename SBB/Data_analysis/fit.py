#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np
from scipy.optimize import leastsq as _leastsq

from SBB.Numpy_extra.numpy_extra import sub_flatten

import pdb

from matplotlib.pyplot import subplots

def polyfit_multi_check(x,Y,Xth=None,deg=1,Y_idxs=None,P=None,ax=None):
    """
    Fast visuallization of a polyfit_multi 
    
    Inputs
    ------
    x     : array like, shape (    M) 
    Y     : array like, shape (...,M) 
    deg   : int
        Degree of the fitting polynomial
    Y_idxs : list of length = ... (aka len(Y_idxs) = Y.ndim-1 )
        The Y idexes to plot
        will be randomly selected by default
    P : The already fitted P 
        will be fitted if not given
    ax : The axis for the plot
    
    Returns
    -------
    ax
    
    See Also
    --------
    polyfit_multi
    
    By default it selects one random index of Y[...,M] (aka a random index in the ... part)
    """
    def get_rdm_idx(n):
        if n == 0 :
            return 0
        else :
            return _np.random.randint(0,n)
    
    if Y_idxs is None :
        # selects a random index to plot
        Y_idxs = [ get_rdm_idx(Y.shape[i]-1) for i in range(Y.ndim -1) ] 
    
    if P is None :
        if Xth is None :
            P = polyfit_multi(x,Y,deg=deg) 
        else :
            P = polyfit_above_th(x,Y,Xth,deg=deg)
    
    if ax is None :
        fig,ax = subplots(1,1)
    
    line, =  ax.plot(x,Y[tuple(Y_idxs)],marker='.',ls='None')
    fit = 0.0
    for i,p in enumerate(P[tuple(Y_idxs)][::-1]) :
        fit += p*x**i
    ax.plot(x,fit,marker='None',ls='--',label='{}'.format(Y_idxs),color=line.get_color())
    ax.legend(title='Y_idxs')
    
    return ax , P

def polyfit_multi(x,Y,deg=1):
    """
    Polynomial fit of Y(X) along the last axis 
    
    Bonus : it works when Y contains nans.

    Inputs
    ------
    X     : array like, shape (    M) 
    Y     : array like, shape (...,M) 
    deg   : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(...,deg+1)
    
    See Also
    --------
    _np.polyfit
    """
    if Y.ndim == 1 :
        Y = Y[None,:]
    
    P_shape = Y.shape[:-1]+(deg+1,)
    P       = _np.full( (_np.prod(P_shape[:-1]),) + P_shape[-1:] , _np.nan )
       
    for j,y in enumerate( sub_flatten(Y) ):  #coulb be made faster by usign sub_flat instead
        not_nan = ~(_np.isnan(y)) 
        all_nan = all( not_nan ==False)
        if all_nan :
            P[j] = _np.nan
        else :
            P[j]   = _np.polyfit(x[not_nan],y[not_nan],deg)
    P.shape = P_shape
    return P

def polyfit_above_th(x,Y,Xth,deg=1):
    """
    For retrocompatibility
    """
    return polyfit_multi_between(x,Y,Xth_low=Xth,Xth_high=None,deg=deg)

def polyfit_multi_between(x,Y,Xth_low=None,Xth_high=None,deg=1):
    """
    Polynomial fit of Y(X) for (X>=Xth_low)&((X=<Xth_high)).

    Inputs
    ------
    X        : array like,                 shape (    M) 
    Y        : array like,                 shape (...,M) 
    Xth_low  : float or array_like with     shape (...  )
    Xth_high : float or array_like with     shape (...  )
    deg      : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(...,deg+1)
    
    See Also
    --------
    _np.polyfit
    """
    if (Xth_low is None) and (Xth_high is None):
        return polyfit_multi(x,Y,deg=deg)
    
    def build_Xth(Xth,Y):
        xth         = _np.full( Y.shape[:-1],_np.nan )
        xth[...]    = Xth
        return xth
        
    def select_x(x,xth_l,xth_h):
        if _np.isnan(xth_l) :
            return x<=xth_h
        elif _np.isnan(xth_h) :
            return x>=xth_l
        else :
            return (x>=xth_l)&(x<=xth_h)
    
    if Y.ndim == 1 :
        Y = Y[None,:]   
    Xth_low  = build_Xth(Xth_low ,Y)  # Memory hungry ... 
    Xth_high = build_Xth(Xth_high,Y)
    
    P_shape = Y.shape[:-1]+(deg+1,)
    P       = _np.full( (_np.prod(P_shape[:-1]),) + P_shape[-1:] , _np.nan )
       
    for j,(y,xth_l,xth_h) in enumerate(zip(sub_flatten(Y),Xth_low.flat,Xth_high.flat)):  #coulb be made faster by usign sub_flat instead
        x_pos = select_x(x,xth_l,xth_h) # Select the abscisse idx repecting the threshold
        not_nan = ~(_np.isnan(y))       
        all_nan = _np.all( (x_pos & not_nan) ==False)
        if all_nan :
            P[j] = _np.nan
        else :
            P[j]   = _np.polyfit(x[x_pos & not_nan],y[x_pos & not_nan],deg)
    P.shape = P_shape
    return P

def lstsq_2D(x0,x1,y,para,fit_func,tol=1e-8,full_output=0): 
    def residuals(p):
        return y-fit_func(x0,x1,p)
    return _leastsq(residuals,para,full_output=full_output,ftol=tol,xtol=tol)
    
def lstsq(x,y,para,fit_func,tol=1e-8,full_output=0): 
    def residuals(p):
        return y-fit_func(x,p)
    return _leastsq(residuals,para,full_output=full_output,ftol=tol,xtol=tol)
    
    