#!/bin/env/python
#! -*- coding: utf-8 -*-
import  numpy as _np

def fourier_transform(F,dt):
    """
    This gives and approximation of F(f) or fourier transform of F(t) if you prefer
    and is not directly the DFT of F(t)    
    """
    return dt*_np.fft.rfft(F)
def ifourier_transform(F,dt,n):
    """
    See Also
    -------
        fourier_transform
    """
    irfft = _np.fft.irfft
    shift = _np.fft.fftshift
    return (1.0/dt)*shift(irfft(F,n=n))

#########################
# Central limit theorem #
#########################

def SE(mu2k,muk,n):
    """ 
        Voir notes Virally Central limit theorem
        Computation of the standard error for the moment of order K
        mu2k : is the moment of order 2 k
        muk  : is the moment of order k
        If these moments are not centered then the definition is good for none centered moment
        Idem for centered moment
    """
    return _np.sqrt(_np.abs(mu2k-muk**2)/float(n))  

#####################
# Moments cumulants #
#####################

def moment(hx,hs,exp,n_total,no_clip=True):
    tmp = ((hx**exp)*hs)
    tmp = tmp[...,1:-1] if no_clip else tmp
    return (tmp.sum(axis=-1))/n_total
   
def centered_moment(hx,hs,exp,n_total,no_clip=True):
    mu = moment(hx,hs,1,n_total,no_clip)
    tmp = (((hx[None,None,:]-mu[...,None])**exp)*hs)
    tmp = tmp[...,1:-1] if no_clip else tmp
    return (tmp.sum(axis=-1))/n_total 

########################
# Numerical derivation #
########################

"""
see: https://web.media.mit.edu/~crtaylor/calculator.html
"""

def compute_differential(X):
    """
        Returns y2-y1 
        in an array of shape (2,shape_of_X_with_a_-1_on_the_last_axis) 
        with [0,...] corresponding to the y2
        and  [1,...] corresponding to the y1
        The user can do y2-y1 afterward to get the differiential
    """
    shape           = X.shape
    shape           = (2,) + shape[:-1] + ( shape[-1]-1, )
    X_diff          = _np.zeros(shape)
    X_diff[0,...]   = X[...,1:]
    X_diff[1,...]   = X[...,:-1]
    return X_diff

###############
# dx variable #
###############

def dydx_centered_3pnts(x,y):
    """
    différence centré pour espacement dx inégale
    
    reference : https://faculty.ksu.edu.sa/sites/default/files/l8_2.pdf
    """
    out = _np.full(y.shape,_np.nan)
    out[...,0]    = (y[...,1]-y[...,0])/(x[...,1]-x[...,0])
    out[...,1:-1] = ( (x[...,1:-1]-x[...,2:])              /((x[...,:-2] -x[...,1:-1])*(x[...,:-2 ]-x[...,2:  ])) )*y[...,0:-2] + \
                    ( (2*x[...,1:-1]-x[...,0:-2]-x[...,2:])/((x[...,1:-1]-x[...,0:-2])*(x[...,1:-1]-x[...,2:  ])) )*y[...,1:-1] +  \
                    ( (x[...,1:-1]-x[...,0:-2])            /((x[...,2:]  -x[...,:-2 ])*(x[...,2:  ]-x[...,1:-1])) )*y[...,2:  ]
    out[...,-1]       = (y[...,-1]-y[...,-2])/(x[...,-1]-x[...,-2])
    return out

##########
# FIX dx #
##########
def central_derivative_3points(dx,y):
    """
    bigO dx**2
    """
    out = _np.zeros(y.shape)
    out[...,1:-1] = ( -y[...,:-2] + y[...,2:]  )/(2.*dx)
    out[...,0] = (y[...,1]-y[...,0])/dx
    out[...,-1]= (y[...,-1]-y[...,-2])/dx
    return out

def central_derivative_9points(dx,y):
    """
    bigO dx**8
    """    
    return (3.*y[:-8]-32.*y[1:-7]+168.*y[2:-6]-672.*y[3:-5]+672.*y[5:-3]-168.*y[6:-2]+32.*y[7:-1]-3.*y[8:])/(840.*dx)
    
 
def central_derivative_21points(dx,y):
    """
    bigO dx**20
    """    
    return (6.192803659266073e+143*y[:-20]-1.3761785840767963e+145*y[1:-19]+1.4707908530326345e+146*y[2:-18]-1.0085422920797508e+147*y[3:-17]+5.00068882203979e+147*y[4:-16]-1.920264487544838e+148*y[5:-15]+6.000826445982782e+148*y[6:-14]-1.6002203612109372e+149*y[7:-13]+3.9005371863629597e+149*y[8:-12]-1.0401434581699147e+150*y[9:-11]-2.3204491725103567e+143*y[10:-10]+1.0401437331570793e+150*y[11:-9]-3.9005386535673158e+149*y[12:-8]+1.6002208363906606e+149*y[13:-7]-6.00082798135261e+148*y[14:-6]+1.9202649089213876e+148*y[15:-5]-5.000689772857675e+147*y[16:-4]+1.0085424604235924e+147*y[17:-3]-1.4707910712381083e+146*y[18:-2]+1.3761787675320434e+145*y[19:-1]-6.192804408049298e+143*y[20:])/(1.1441580589981547e+150*dx)
 
 
def central_derivative_kernel(dx,bigO=2):
    """
    This can be automated to any order of derivative and bigO.
    https://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0935077-0/S0025-5718-1988-0935077-0.pdf
    """
    if bigO==2  :
        return r_[-1.,0,1]/(2.*dx)
    elif bigO==8 :
        return r_[3.,-32.,168.,-672.,0,672.,-168.,32.,-3.]/(840.*dx)
    elif bigO==20 :
        return r_[6.192803659266073e+143,-1.3761785840767963e+145,1.4707908530326345e+146,-1.0085422920797508e+147,5.00068882203979e+147,-1.920264487544838e+148,6.000826445982782e+148,-1.6002203612109372e+149,3.9005371863629597e+149,-1.0401434581699147e+150,-2.3204491725103567e+143,1.0401437331570793e+150,-3.9005386535673158e+149,1.6002208363906606e+149,-6.00082798135261e+148,1.9202649089213876e+148,-5.000689772857675e+147,1.0085424604235924e+147,-1.4707910712381083e+146,1.3761787675320434e+145,-6.192804408049298e+143]/(1.1441580589981547e+150*dx)
    else:
        return r_[0.,]