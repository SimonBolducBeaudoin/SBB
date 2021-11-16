#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy
from SBB.Utilities.General_tools import symetrize,fourier_transform,find_nearest_A_to_a
from SBB.Utilities import Noise_Theory_Junction 

def binV2_to_A2(S2,R_acq,mv_per_bin):
    """
    Used to convert SII(t)[bin_V**2] to SII(t)[A**2]
    """
    return S2*(mv_per_bin*1.0e-3)**2/(R_acq**2)
    
def window_after_2ns(S2):
        """
            Damping everything more than 2 ns
            At a sampling rate of 0.03125 it means everything after the 64th point
            
            This will need to be rewritten if used with another aquisition card...
        """ 
        def damp(x,epsilon,x_0):
            return numpy.exp((-1)*epsilon*(x-x_0))
        def compute_epsilon(red,after_lenght):
            return -numpy.log(1.0/red)/(after_lenght)
        red = 1000
        L_0 = 65
        epsilon = compute_epsilon(red,after_lenght=L_0)
        shape = S2.shape
        len = shape[-1]
        out = numpy.zeros(shape)
        out = S2
        for index in range(L_0,len):
            out[...,index] = S2[...,index]*damp(index,epsilon,L_0-1)
        return out

def SII_dc_of_t_to_spectrum(S2,dt):
    S2_windowed       = window_after_2ns(S2)
    S2_sym            = symetrize (S2_windowed)
    return numpy.abs(fourier_transform(S2_sym,dt))
    
def compute_Ith(f,f_max,Te,eps=0.01,R_jct=50.0):
    """
    Computes the threshold currents above which SII = eI in good approximation for each frequencies.

    This function intended to help with the use of compute_G_and_SII_amp

    Parameters
    ----------
    f[Hz]       : float or array
    f_max[Hz]   : float
        frequencies above this value will get the same threshold as this frequency
    Te [K]      : float
        expected temperature of the electron of the junction
    eps [~]     : float
        Error on the convergence of Seq
    R_jct [Ohm] : float
    """
    V_th = Noise_Theory_Junction.V_th(f,Te=Te,epsilon=eps)
    _,f_max_idx = find_nearest_A_to_a(f_max,f)
    f_max_idx = f_max_idx[0]
    V_th[f_max_idx:] = V_th[f_max_idx]
    return V_th/R_jct

def compute_noiseTemp(SII,G,R=50.0):
    """
    Used to compute the amplifier's noise temperature
    Converts SII[A**2/Hz] in K

    Last update
    -----------
    Created
    I choose the factor 2.0 instead of 4.0 in Jonhson-Nyquist Noise
    """
    return SII*R/(2.0*C.k*G)

 