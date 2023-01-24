#!/bin/env/python
#! -*- coding: utf-8 -*-

from numpy import exp,vectorize,linspace,zeros, r_
from SBB.Utilities.Memoize import MemoizeMutable
from scipy.special import iv,comb # modified bessels
from scipy.optimize import brentq

def _gaussian_2D(x,y,sigma_x=1.,sigma_y=1.,rho=0.):
    "General 2d-gaussian"
    prefac = 1./(2.*pi*sigma_x*sigma_y*sqrt(1.-rho**2))
    return prefac * exp( (-1./(2*(1-rho**2))) * ( (x/sigma_x)**2 - 2.*rho*(x/sigma_x)*(y/sigma_y) + (y/sigma_y)**2 )  )
gaussian_2D = MemoizeMutable(vectorize(_gaussian_2D))


def _rotating_gaussian(r, R):
    """
    Not normalized !!!
    e^{-r^2} I_0 \bigg( r^2\frac{1-R^2}{1+R^2} \bigg)
    """
    A = r**2*(1.-R**2)/(1.+R**2)
    if A > 713 :
        """
        I assume that the exp dominates 
        """
        return 0.
    else :
        return exp( -r**2 ) * iv(0, A)
    
rotating_gaussian = MemoizeMutable(vectorize(_rotating_gaussian))
    

def moment_of_f(r,f_of_r,power) :
    return (r**power * f_of_r ).sum()
    
def std_moment_of_f(r,f_of_r,power) :
    if power == 1 :
        return 0
    elif power == 2:
        return 1
    else :
        var  = moment_of_f(r,f_of_r,2).sum()
        mu_k = moment_of_f(r,f_of_r,power).sum()
        return mu_k/var**(power/2.0)
    
def std_cumulants_recursive(r,f_of_r,n):
    if n == 1 :
        return 0 
    elif n == 2 :
        return 1
    else :
        m_n = std_moment_of_f(r,f_of_r,n)
        sum = 0.0
        for k in range(1,n) :
            sum += comb(n-1,k-1)*std_cumulants_recursive(r,f_of_r,k) * std_moment_of_f(r,f_of_r,n-k)
        return m_n - sum

def C4_of_R(R,r_end=50.,r_len = 10000):
    """
    Invalid when R ==> 0
    """
    r = linspace(0,r_end,r_len)
    r  = r_[-1.*r[::-1],r] #symmetrized 
    f_of_r    = rotating_gaussian(r,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    f_of_r    = f_of_r/norm
    r2        = moment_of_f(r,f_of_r,2)
    r4        = moment_of_f(r,f_of_r,4)
    return r4-3.*r2**2

def C2_of_R(R,r_end=50.,r_len = 100000):
    """
    Invalid when R ==> 0
    """
    r = linspace(0,r_end,r_len)
    r  = r_[-1.*r[::-1],r] #symmetrized 
    f_of_r    = rotating_gaussian(r,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    f_of_r    = f_of_r/norm
    return  moment_of_f(r,f_of_r,2)

def C4_over_C2_2_of_R(R,r_end=50.,r_len = 10000):
    """
    Invalid when R ==> 0
    """
    r = linspace(0,r_end,r_len)
    r  = r_[-1.*r[::-1],r] #symmetrized 
    f_of_r    = rotating_gaussian(r,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    f_of_r    = f_of_r/norm
    r2        = moment_of_f(r,f_of_r,2)
    r4        = moment_of_f(r,f_of_r,4)
    return r4/r2**2-3.
    

def zero_C4_over_C2_2(R,y,r_end=50.0,r_len = 10000):
    return C4_over_C2_2_of_R(R,r_end,r_len) - y

def find_R_from_C4_over_C2_2(C4_norms,maxiter=100, xtol= 1.e-3, rtol=1.e-3,r_end=50.0,r_len = 10000):
    R = zeros(C4_norms.shape)
    
    try :
        def factory(R,y,r_end=r_end,r_len = r_len):
            return C4_over_C2_2_of_R(R,r_end,r_len) - y   
        for i,y in enumerate (C4_norms) :
            if y <= 0.01:   # for convenience
                R[i] = 1.0
            else :
                R[i] = brentq(factory, 0.01, 1., args=(y,) ,maxiter=maxiter, xtol=xtol, rtol=rtol,full_output=False, disp=True)
    except ValueError:
        def factory(R,y,r_end=r_end,r_len = r_len*10): # Augmenting the number of points if brentq show an error for the initial args
            return C4_over_C2_2_of_R(R,r_end,r_len) - y   
        for i,y in enumerate (C4_norms) :
            if y <= 0.01:   # for convenience
                R[i] = 1.0
            else :
                R[i] = brentq(factory, 0.01, 1., args=(y,) ,maxiter=maxiter, xtol=xtol, rtol=rtol,full_output=False, disp=True)
    return R


def C4_over_C2_of_R(R,r_end=50.,r_len = 10000):
    """
    Computes (<r**4> - 3 <r**2>) /<r**2>
    """
    r = linspace(0,r_end,r_len)
    r  = r_[-1.*r[::-1],r] #symmetrized 
    f_of_r    = rotating_gaussian(r,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    f_of_r    = f_of_r/norm
    r2        = moment_of_f(r,f_of_r,2)
    r4        = moment_of_f(r,f_of_r,4)
    return (r4-3.*r2**2)/r2
    
def C2_over_C0_of_R(R,r_end=50.,r_len = 10000):
    """
    Computes ( <r**2>) /<r**0>)
    """
    r = linspace(0,r_end,r_len)
    r  = r_[-1.*r[::-1],r] #symmetrized 
    f_of_r    = rotating_gaussian(r,R=R)
    norm      = moment_of_f(r,f_of_r,0)
    r2        = moment_of_f(r,f_of_r,2)
    return r2/norm
 
# def find_sigma_square(C4,C2,R,r_end=50.,r_len = 10000):
    # tmp = zeros(R.shape)
    # for i,RR in enumerate(R) :
        # tmp[i] = C4_over_C2_of_R(RR,r_end=50.,r_len = 10000)
    # return 0.25* C4 / (C2 * tmp) 
 
def find_sigma_square(C2,R,r_end=50.,r_len = 10000):
    tmp = zeros(R.shape)
    for i,RR in enumerate(R) :
        tmp[i] = C2_over_C0_of_R(RR,r_end=50.,r_len = 10000)
    return 0.25* C2 / tmp 
    