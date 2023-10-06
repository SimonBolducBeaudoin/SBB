#!/bin/env/python
#! -*- coding: utf-8 -*-

"""
    This module was written to modelize rotating_gaussians as one would expect form a rotating squeezed distribution.
"""

import numpy as _np
import scipy.integrate as _int
from scipy.special import iv as _iv # modified bessels
from scipy.optimize import brentq as _brentq

def _rotating_gauss(r, R):
    """
    Not normalized !!!
    e^{-r^2} I_0 \bigg( r^2\frac{1-R^2}{1+R^2} \bigg)
    """
    A = r**2*(1.-R**2)/(1.+R**2)
    if A > 713:
        """
        I assume that the exp dominates 
        """
        return 0.
    elif r==0 :
        return 1.
    else :
        return _np.exp( -r**2 ) * _iv(0, A)
rotating_gauss = _np.vectorize(_rotating_gauss)

def _expval_r2k(k,R,r_end=_np.inf):
    """
    returns <r^2k>/<r^0>
    """
    r0  =  _int.quad( lambda r: _rotating_gauss(r,R) , 0, r_end,full_output=True)
    r2k =  _int.quad( lambda r: r**2*_rotating_gauss(r,R) , 0, r_end,full_output=True)
    return r2k[0]/r0[0]
expval_r2k = _np.vectorize(_expval_r2k,excluded=['r_end'])

def _get_varx_vary(C2,R,r_end=_np.inf):
    """
    returns sigma_x^2 and sigma_y^2 
    """
    r2 = _expval_r2k(1,R,r_end)
    return 0.25*C2*(1.0+R**2)/r2 ,0.25*C2*(1.0/R**2+1.0)/r2
get_varx_vary = _np.vectorize(_get_varx_vary,excluded=['r_end'])
    
def _Std_C4(R,r_end=_np.inf):
    r0      =  _int.quad( lambda r: _rotating_gauss(r,R) , 0, r_end,full_output=True)
    r2      =  _int.quad( lambda r: r**2*_rotating_gauss(r,R) , 0, r_end,full_output=True)
    r4      =  _int.quad( lambda r: r**4*_rotating_gauss(r,R) , 0, r_end,full_output=True)
    return r0[0]*r4[0]/r2[0]**2-3.
Std_C4 = _np.vectorize(_Std_C4,excluded=['r_end'])
    
def _find_R(StdC4,maxiter=100, xtol= 1.e-3, rtol=1.e-3,r_end=_np.inf,force_physical_C4=True):
    """
    Tries to find R given StdC4 using Brent’s method (root finding method).
    """
    if force_physical_C4 :
        if StdC4< 0  :
            return 1.0
    def residual(R,y,r_end=r_end):
        return _Std_C4(R,r_end) - y 
    try :
        return _brentq(residual, 0.0, 1., args=(StdC4,) ,maxiter=maxiter, xtol=xtol, rtol=rtol,full_output=False, disp=True)
    except :
        return _np.nan
find_R = _np.vectorize(_find_R,excluded=['maxiter','xtol','rtol','r_end'])

    