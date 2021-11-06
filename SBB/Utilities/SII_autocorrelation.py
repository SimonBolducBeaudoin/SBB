#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy

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
            out[:,index] = S2[:,index]*damp(index,epsilon,L_0-1)
        return out
        
 