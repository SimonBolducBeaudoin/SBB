#!/bin/env/python
#! -*- coding: utf-8 -*-

from numpy import exp,vectorize,linspace,zeros
from Memoize import MemoizeMutable
from scipy.special import iv # modified bessels
from scipy.optimize import brentq

def _gaussian_2D(x,y,sigma_x=1.,sigma_y=1.,rho=0.):
    "General 2d-gaussian"
    prefac = 1./(2.*pi*sigma_x*sigma_y*sqrt(1.-rho**2))
    return prefac * exp( (-1./(2*(1-rho**2))) * ( (x/sigma_x)**2 - 2.*rho*(x/sigma_x)*(y/sigma_y) + (y/sigma_y)**2 )  )
gaussian_2D = MemoizeMutable(vectorize(_gaussian_2D))

def _rotating_gaussian(r,sigma_x = 1., R = 1.):
    return 0.5* exp( -0.25*r**2*( 1 + R**2 )/sigma_x**2 ) * iv(0, 0.25*r**2*( 1 - R**2 )/sigma_x**2 )
rotating_gaussian = MemoizeMutable(vectorize(_rotating_gaussian))

def moment_of_f(r,f_of_r,power) :
    return (r**power * f_of_r ).sum()
    
def C2_norm_vs_R(R,sigma_x=1.):
    """
    Numerically invalid when R ==> 0
    C2/C2(R=1)
    Note that sigma_x**2 is the variance of the gaussian case R=1.
        SO sigma_x**2 = C2(R=1)
    """
    r = linspace(0,50*sigma_x,10000)
    f_of_r    = rotating_gaussian(r,sigma_x=sigma_x,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    r2        = moment_of_f(r,f_of_r,2)/norm
    return r2/sigma_x**2

def zero_C2_over_C2ref(R,y):
    return C2_norm_vs_R(R) - y

def find_R_from_C2_over_C2ref(C2_over_C2ref):
    R = zeros(C2_over_C2ref.shape)
    for i,y in enumerate (C2_over_C2ref) :
        R[i] = brentq(zero_C2_over_C2ref, 0.1, 1., args=(y,) ,maxiter=100, xtol= 1.e-3, rtol=1.e-3,full_output=False, disp=True)
    return R

def C4_norm_vs_R(R,sigma_x=1.):
    """
    Numerically invalid when R ==> 0
    mu4/mu2**2 - 3
    """
    if R==1:
        return 0.
    r = linspace(0,50*sigma_x,10000)
    f_of_r    = rotating_gaussian(r,sigma_x=sigma_x,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    r2        = moment_of_f(r,f_of_r,2)/norm
    r4        = moment_of_f(r,f_of_r,4)/norm
    return r4/(r2**2) - 3

def zero_C4_norms(R,y):
    return C4_norm_vs_R(R) - y

def find_R_from_C4_norms(C4_norms,maxiter=100, xtol= 1.e-3, rtol=1.e-3):
    R = zeros(C4_norms.shape)
    for i,y in enumerate (C4_norms) :
        if y <= 0:
            R[i] = 1.0
        else :
            R[i] = brentq(zero_C4_norms, 0.1, 1., args=(y,) ,maxiter=maxiter, xtol=xtol, rtol=rtol,full_output=False, disp=True)
    return R
