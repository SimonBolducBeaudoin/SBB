#!/bin/env/python
#! -*- coding: utf-8 -*-

from scipy import tanh as _tanh , arctanh as _arctanh
import numpy as _np
from scipy.special import jv as _besselJ
import scipy.constants as _C

from SBB.Brrrr.Memoize import MemoizeMutable as _MemMut

###############
# Description #
###############

"""
Theoretical autocovariance/noise fonctions.

- We use analitical limits when possible for special cases.
    - We should only get ±inf when it's the actual analitical result.
- We define everything in its proper representation; no FFTs are used.
- Vectorization is used for convenience with special cases.
    - Code should still be performant enough for live plotting.
    - Memoization makes adjusting graphics much faster.

Default units are A² for autocovariance and A²/Hz for spectral densities.
    - R can be set to 2*C.k to obtain Kelvin.
    - Or you can multiply the A²/Hz result by R/(2*C.k).
"""


#################################################################################
def coth(x):
    return 1./_tanh(x)

# Optimizing this function helps a lot in the end.
def xcothx(x):
    x = _np.r_[x]
    ret = _np.empty(x.shape)
    mask = x==0
    xmask = x[~mask]
    ret[mask]=1
    ret[~mask]=xmask*coth(xmask)
    return ret

# Equivalent and clearer but slower implementation of xcothx
@_np.vectorize
def _xcothx(x):
    if x==0:
        return 1
    else:
        return x*coth(x)

################################################################

def Vsquare_to_K(SII,Z_jct):
    """ 
    Converts from V**2/Hz to Kelvin
    Si on replie les fréquences negatives sur les fréquences positives 
    on prend le facteur 1/2K (i.e. l'abscisse n'existe pas sur les fréquence negative)
    Sinon on prend 1/4K
    """
    return SII/(2.0*Z_jct*_kb)

def K_to_Vsquare(SII,Z_jct):
    """ 
    Converts from Kelvin to V**2/Hz
    """
    return (2.0*Z_jct*_kb)*SII
    
def Adimensional_to_K(SII,Te):
    """
    Converts a SII with no unit to Kelvin
    """
    return Te*SII

def Adimensional_to_Vsquare(SII,Z_jct,Te):
    """
    Converts a SII with no unit to V**2/Hz
    """
    return (2.0*Z_jct*_kb*Te)*SII

def Adimensional_to_Asquare(SII,Z_jct,Te):
    """
    Converts a SII with no unit to A**2/Hz
    """
    return (2.0*_kb*Te)*SII/Z_jct

def V_th(f,Te=None,epsilon=0.01):
    """
        Return the V_th for which e*V_th =  h f
        
        f : float or 1D array
        
        When Te != None
            Returns a more conservative value for V_th
            for which SII_eq = asymtotic value + epsilon
    """
    if Te :
        cst = _arctanh(1.0/(1.0+epsilon)) # coth(cst) = 1.01
        return _np.max(_np.array([  _np.abs( cst*2.0*_C.k *Te - _C.h*f ) ,  _np.abs( cst*2.0*_C.k*Te + _C.h*f ) ]),axis=0)/_C.e
    else :
        return _C.h*f/_C.e

#####################################################################

def compute_nu(V):
    return _C.e*V/_C.hbar

def Seq_of_t(tau,Te,R):
    """
    Python implementation of Seq(tau)
    
    Time domain auto-correlation of a tunnel junction at thermal equilibrium
    
    Latex :
        S_{eq} = - \frac{\pi (k T)^2}{R \hbar} \_np.sinh^{-2} \bigg( \frac{\pi k T}{\hbar} \tau \bigg)
        
    See Also
    --------
    Thesis : Mesures temporelles large bande ..., Simoneau , equation # 7.57
        
    Parameters
    ----------
    tau : float     Time delay [s]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2]       
    """
    if Te==0:
        return -1./tau**2*_C.hbar/(pi*R)
    else:
        return -pi*(_C.k*Te)**2/(R*_C.hbar)*1./(_np.sinh(pi*_C.k*Te*tau/_C.hbar))**2
        
def _Sdc_of_t(tau,nu,Te,R):
    """
    Time domain auto-correlation of a tunnel junction at Temperature T and DC bias Vdc
    
    Parameters
    ----------
    tau : float     Time delay [s]
    nu  : float     eV/hbar [s-1]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2]
    """
    return Seq_of_t(tau,Te,R)*_np.cos(nu*tau)
Sdc_of_t = _MemMut(_np.vectorize(_Sdc_of_t))

def _DSdc_of_t(tau,nu,Te,R):
    """
    Contribution to the auto-correlation of a tunnel junction at thermal equilibrium due to DC bias 
    Sometime called excess noise
    
    Parameters
    ----------
    tau : float     Time delay [s]
    nu  : float     eV/hbar [s-1]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2]
    """
    return -2*Seq_of_t(tau,Te,R)*_np.sin(nu*tau/2.)**2
    # Equivalent form for testing
    #return _tSdc(tau,nu,Te,R)-_tSdc(tau,0,Te,R)
DSdc_of_t = _MemMut(_np.vectorize(_DSdc_of_t))
 
 
def _SPH_of_t(tau,nu,Te,R):
    """
    Temperature excess noise or Pauli-Heisenberg oscillations
    
    Parameters
    ----------
    tau : float     Time delay [s]
    nu  : float     eV/hbar [s-1]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2]
    """
    if tau==0:
        D = pi*(_C.k*Te)**2/(3*_C.hbar*R)
    else:
        D = Seq_of_t(tau,Te,R)-Seq_of_t(tau,0,R)
    return D*_np.cos(nu*tau)
SPH_of_t = _MemMut(_np.vectorize(_SPH_of_t))


####################
# Frequency domain #
####################

def _Seq_of_f(omega,Te,R):
    """
    Python implementation of Seq(tau)
    
    Frequency domain auto-correlation of a tunnel junction at thermal equilibrium
    
    Latex :
        S_{eq}(\omega) = \frac{2kT}{R} \frac{\hbar\omega}{2kT} \coth \frac{\hbar \omega}{ 2 k T }
        
    See Also
    --------
    Thesis : Mesures temporelles large bande ..., Simoneau , equation # 7.48
        
    Parameters
    ----------
    tau : float     Time delay [s]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2/Hz]       
    """
    if Te==0:
        return abs(_C.hbar*omega)/R
    else:
        return 2*_C.k*Te/R * xcothx(_C.hbar*omega/(2*_C.k*Te))
Seq_of_f = _MemMut(_np.vectorize(_Seq_of_f))

def _Sdc_of_f(omega,nu,Te,R):
    return (_Seq_of_f(nu-omega,Te,R)+_Seq_of_f(nu+omega,Te,R))/2.
Sdc_of_f = _MemMut(_np.vectorize(_Sdc_of_f))

def bessel_weights(n,z):
    return _besselJ(Ns,z)**2.
    
def _Spa_of_f(omega,nu,nuac,Omega,Te,R,nBessel=21):
    if not nuac*Omega:
        return _Sdc_of_f(omega,nu,Te,R)
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = arange(-nBessel//2+1,nBessel//2+1)
    freqs = omega+Ns*Omega

    Sdcs = _Sdc_of_f(freqs,nu,Te,R)
    bessels = _besselJ(Ns,z)**2.

    return _np.dot(bessels, Sdcs)
Spa_of_f = _MemMut(_np.vectorize(_Spa_of_f))

def _Sphi_of_f(phi,omega,nu,nuac,Omega,Te,R,nBessel=21):
    if phi is None or not nuac*Omega:
        return _Spa(omega,nu,Te,nuac,Omega,R,nBessel=nBessel)
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = arange(-nBessel//2+1,nBessel//2+1)
    freqs_m = -omega+Ns*Omega+nu
    freqs_p = +omega+Ns*Omega+nu

    Sds_m = _Seq_of_f(freqs_m,Te,R)*exp(+1j*Ns*phi)
    Sds_p = _Seq_of_f(freqs_p,Te,R)*exp(-1j*Ns*phi)

    bessels = _besselJ(Ns,z)

    return 0.5*(exp(-1j*z*_np.sin(phi))*_np.dot(bessels,Sds_m) + exp(+1j*z*_np.sin(phi))*_np.dot(bessels,Sds_p))
Sphi_of_f = _MemMut(_np.vectorize(_Sphi_of_f))

def _Spa_adiabatique(omega, phi, nu, nuac, Te, R) :
    dphi  = phi[1] - phi[0]
    return (1./(2*pi)) * Sdc_of_f( omega, nu + nuac , Te , R ).sum(axis=0) * dphi
_Spa_adiabatique = _MemMut(_np.vectorize(_Spa_adiabatique))

def _betap(p,f,nu,Te,nuac,Omega,R,nBessel=21):
    omega = 2.0*pi*f
    if p==0:
        return _Spa(omega,nu,Te,nuac,Omega,R,nBessel=nBessel)
    if not nuac*Omega:
        return 0
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = arange(-nBessel//2+1,nBessel//2+1)
    freqs_m = -omega-Ns*Omega+nu
    freqs_p = +omega+Ns*Omega+nu

    Sds_m = _Seq(freqs_m,Te,R)*(-1.)**p
    Sds_p = _Seq(freqs_p,Te,R)

    bessels = _besselJ(Ns,z)*_besselJ(Ns+p,z)

    return _np.dot(bessels, (Sds_m+Sds_p)/2.)
betap = _MemMut(_np.vectorize(_betap))

def Xp(p,omega,nu,Te,nuac,Omega,R,nBessel=21):
    return betap(p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)+betap(-p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)

def Yp(p,omega,nu,Te,nuac,Omega,R,nBessel=21):
    return betap(p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)-betap(-p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)

def Sqz(p,omega,nu,Te,nuac,Omega, R,nBessel=21):
    return _np.sign(p)*betap(abs(p),omega,nu,Te,nuac,Omega,R,nBessel=21)+Spa(omega,nu,Te,nuac,Omega,R,nBessel=21)