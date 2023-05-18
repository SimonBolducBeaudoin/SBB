#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np

def gamma_to_Z(gamma,Z0=50.0):
    """
        Converts reflexion coefficient to impedance
        gamma = ZL - Z0 /  ZL + Z0*
    """
    return Z0*(gamma+1.0)/(1.0-gamma)
    
def Z_to_gamma(Z,Z0=50.0):
    return (Z-Z0)/(Z+Z0)
    
def dB_to_V2_over_V1(dB,R2,R1=50.0):
    """
    V2/V1  = 10**(dBm/20) sqrt(R2/R1)
    """
    return _np.sqrt(R2/R1)*10.0**(dB/20.0)
    
def dBm_to_V(dBm,R2,R1=50.0):
    """
        Returns the peak voltage from dBm
        
        dBm  = 10 log10( V2**2  R1 )
                         ------ --
                         V1**2  R2
        V2  = 10**(dBm/20) sqrt(2 R2 V1 )
        sqrt(2) converts rms to peak, 
        R1 = 50.0 (in general), 
        V1 = 0.001[W]
    """    
    return _np.sqrt(2.0*R2*0.001)*10.0**(dBm/20.0)

def lin_to_dB(x):
    return 10*_np.log10(x)

def dB_to_lin(x):
    return 10.0**(x/10.0)