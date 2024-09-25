#!/bin/env/python
#! -*- coding: utf-8 -*-

from __future__ import division
from past.utils import old_div
from .Tunnel_Junction import *

######################
# Photon experiments #
######################

def C2_of_f(omega,nu,nuac=0.,Omega=0.,Te=0.050,R=50.0,nBessel=21):
    """
    Voir : Simoneau et al. Photon-pair shot noise in electron shot noise
    """
    if (nuac == 0.) or (Omega == 0.):
        return (old_div(R,(_C.hbar*omega)))*Sdc_of_f(omega,nu,Te,R) * 0.5
    else :
        z = old_div(nuac,Omega)
        nBessel -= nBessel%2-1    # Ensure it's odd
        Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
        Sum = 0.
        for n in Ns:
            Sum +=  0.5*_besselJ(n,z)*_besselJ(n,z)*( Seq_of_f(nu+omega+n*Omega,Te,R) + Seq_of_f(nu-omega-n*Omega,Te,R) ) 
        return (old_div(R,(_C.hbar*omega)))*Sum.sum() * 0.5 # this is necessary 
        
def dC2ac_th(omega,nu,nuac=0.,Omega=0.,Te=0.050,R=50.0,nBessel=21):
    return C2_of_f(omega,nu,nuac,Omega,Te,R,nBessel) - C2_of_f(omega,nu,0.,Omega,Te,R,nBessel)

def n_th(omega,nu,nuac=0.,Omega=1.e9,Te=0.050,R=50.0,nBessel=21):
    return C2_of_f(omega,nu,nuac,Omega,Te,R,nBessel) - 0.5 

def dnac_th(omega,nu,nuac=0.0,Omega=1e9,Te=0.050,R=50.0,nBessel=21):
    return n_th(omega,nu,nuac,Omega,Te,R,nBessel) - n_th(omega,nu,0.,Omega,Te,R,nBessel)

def C4_th(omega,nu,nuac,Omega,Te=0.050,R=50.0,nBessel=21):
    """
    Voir : Simoneau et al. Photon-pair shot noise in electron shot noise
    """
    if omega == old_div(Omega,2) : # dangerous
        z = old_div(nuac,Omega)
        nBessel -= nBessel%2-1    # Ensure it's odd
        Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
        Sum = 0.
        for n in Ns:
            Sum +=  0.5*_besselJ(n,z)*_besselJ(n+1,z)*( Seq_of_f(nu+omega+n*Omega,Te,R) - Seq_of_f(omega-n*Omega-nu,Te,R) ) 
        return ((3./2.)*((old_div(R,(_C.hbar*omega)))*Sum.sum())**2) * 0.25 # this is necessary 
    else :
        return 0.
#C4_th = _MemMut(_np.vectorize(_C4_th))

#def C4ac_th(omega,nu,nuac,Omega,Te=Te,R=R,nBessel=21):
#    return C4_th(omega,nu,nuac,Omega,Te,R,nBessel) - C4_th(omega,nu,0.,Omega,Te,R,nBessel)

def dn2_th(omega,nu,nuac=None,Omega=1.e9,Te=0.050,R=50.0,nBessel=21):
    """
    dn2 = 2/3 C_4 + C_2^2 - 1/4
    """
    if ( (nuac == None) or (nuac == 0.0) )   :
        return C2_th(omega,nu,0.0,Omega,Te,R,nBessel)**2 - 0.25 # this is an attempt
    else :
        return 2./3. * C4_th(omega,nu,nuac,Omega,Te,R,nBessel) + C2_th(omega,nu,nuac,Omega,Te,R,nBessel)**2 - 0.25

def ddn2ac_th(omega,nu,nuac=0.,Omega=1.e9,Te=0.050,R=50.0,nBessel=21):
    return 2./3. * C4_th(omega,nu,nuac,Omega,Te,R,nBessel) + nac_th(omega,nu,nuac,Omega,Te,R,nBessel)*(nac_th(omega,nu,nuac,Omega,Te,R,nBessel)+1)

def Fano_th(omega,nu,nuac,Omega,Te=0.050,R=50.0,nBessel=21):
    return old_div(dn2_th(omega,nu,nuac,Omega,Te,R,nBessel),n_th(omega,nu,nuac,Omega,Te,R,nBessel))

def Fanoac_th(omega,nu,nuac,Omega,Te=0.050,R=50.0,nBessel=21):
    return old_div(dn2ac_th(omega,nu,nuac,Omega,Te,R,nBessel),nac_th(omega,nu,nuac,Omega,Te,R,nBessel))