#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np
from scipy.special import jv as _besselJ
import scipy.constants as _C

from SBB.Brrrr.Memoize import MemoizeMutable as _MemMut

from numba import njit,vectorize, float64 , int32, int64

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

@vectorize([float64(float64,float64,float64)])
def R_coulomb_block(I,Vo,D):
    """
    A resistance model when there is coulomb blockade
    V = RI + Vo tanh(I/D)
    Return Vo tanh(I/D)/I ==> The difference in resistance due to the coulomb blockade
    tanh(x)/x when x = 0 ==> 1 
    """
    if I == 0.0:
        return Vo/abs(D)
    else :
        return Vo*_np.tanh(I/abs(D))/I
        
#################################################################################
@vectorize([float64(float64),float64(int64)])
def coth(x):
    return 1./_np.tanh(x)

# Optimizing this function helps a lot in the end.
@vectorize([float64(float64),float64(int64)])
def xcothx(x):
    return x/_np.tanh(x) if x!=0 else 1.0

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
        cst = _np.arctanh(1.0/(1.0+epsilon)) # coth(cst) = 1.01
        return _np.max(_np.array([  _np.abs( cst*2.0*_C.k *Te - _C.h*f ) ,  _np.abs( cst*2.0*_C.k*Te + _C.h*f ) ]),axis=0)/_C.e
    else :
        return _C.h*f/_C.e

#####################################################################

def compute_nu(V):
    return _C.e*V/_C.hbar
    
@vectorize([float64(float64,float64,float64)])
def Seq_of_t(tau,Te,R):
    """
    Seq_of_t(tau,Te,R)
    Time domain auto-correlation of a tunnel junction at thermal equilibrium
    S_{eq} = - \frac{\pi (k T)^2}{R \hbar} \_np.sinh^{-2} \bigg( \frac{\pi k T}{\hbar} \tau \bigg)
        
    Parameters
    ----------
    tau : float     Time delay [s]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2]       
    
    See Also
    --------
    Thesis : Mesures temporelles large bande ..., Simoneau , equation # 7.57
    
    """
    if Te==0:
        return -1./tau**2*_C.hbar/(_np.pi*R)
    else:
        return -_np.pi*(_C.k*Te)**2/(R*_C.hbar)*1./(_np.sinh(_np.pi*_C.k*Te*tau/_C.hbar))**2

@vectorize([float64(float64,float64,float64,float64)])    
def Sdc_of_t(tau,nu,Te,R):
    """
    Sdc_of_t(tau,nu,Te,R)
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
#Sdc_of_t = _MemMut(_np.vectorize(_Sdc_of_t))

@vectorize([float64(float64,float64,float64,float64)])    
def DSdc_of_t(tau,nu,Te,R):
    """
    DSdc_of_t(tau,nu,Te,R)
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
#DSdc_of_t = _MemMut(_np.vectorize(_DSdc_of_t))
 
@vectorize([float64(float64,float64,float64,float64)])    
def SPH_of_t(tau,nu,Te,R):
    """
    SPH_of_t(tau,nu,Te,R)
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
        D = _np.pi*(_C.k*Te)**2/(3*_C.hbar*R)
    else:
        D = Seq_of_t(tau,Te,R)-Seq_of_t(tau,0,R)
    return D*_np.cos(nu*tau)
#SPH_of_t = _MemMut(_np.vectorize(_SPH_of_t))

####################
# Frequency domain #
####################

@vectorize([float64(float64,float64,float64)])   
def Seq_of_f(omega,Te,R):
    """
    Seq_of_f(omega,Te,R)
    Frequency domain auto-correlation of a tunnel junction at thermal equilibrium
    S_{eq}(\omega) = \frac{2kT}{R} \frac{\hbar\omega}{2kT} \coth \frac{\hbar \omega}{ 2 k T }    
    
    Parameters
    ----------
    omega : float     Time delay [s]
    Te  : float     Electron Temperature [K]
    R   : float     Junction resistance [Ohm]
    Returns
    -------
    Seq : float     [A^2/Hz]   

    See Also
    --------
    Thesis : Mesures temporelles large bande ..., Simoneau , equation # 7.48    
    """
    if Te==0:
        return abs(_C.hbar*omega)/R
    else:
        return 2*_C.k*Te/R * xcothx(_C.hbar*omega/(2*_C.k*Te))
#Seq_of_f = _MemMut(_np.vectorize(_Seq_of_f))

@vectorize([float64(float64,float64,float64,float64)])   
def Sdc_of_f(omega,nu,Te,R):
    """
    Sdc_of_f(omega,nu,Te,R)
    """
    return (Seq_of_f(nu-omega,Te,R)+Seq_of_f(nu+omega,Te,R))/2.
#Sdc_of_f = _MemMut(_np.vectorize(_Sdc_of_f))

@vectorize([float64(float64,float64,float64,float64,float64)])   
def Sdc_asym_of_f(omega,nu,nu0,Te,R):
    """
    Asymetric noise due to voltage offset ( ex :  if there is presence of thermoelectric effect )
    Sdc_of_f(omega,nu,nu0,Te,R)
    i.e. ( S(V-V0) - S(V+V0) )/2
    """
    return ( Sdc_of_f(omega,nu-nu0,Te,R)-Sdc_of_f(omega,nu+nu0,Te,R) )/2.0

def bessel_weights(n,z):
    return _besselJ(Ns,z)**2.

def _Spa_of_f(omega,nu,nuac,Omega,Te,R,nBessel=21):
    if not nuac*Omega:
        return Sdc_of_f(omega,nu,Te,R)
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
    freqs = omega+Ns*Omega

    Sdcs = Sdc_of_f(freqs,nu,Te,R)
    bessels = _besselJ(Ns,z)**2.
    return _np.dot(bessels, Sdcs)
Spa_of_f = _MemMut(_np.vectorize(_Spa_of_f))

def _Sphi_of_f(phi,omega,nu,nuac,Omega,Te,R,nBessel=21):
    if phi is None or not nuac*Omega:
        return Spa(omega,nu,Te,nuac,Omega,R,nBessel=nBessel)
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
    freqs_m = -omega+Ns*Omega+nu
    freqs_p = +omega+Ns*Omega+nu

    Sds_m = Seq_of_f(freqs_m,Te,R)*exp(+1j*Ns*phi)
    Sds_p = Seq_of_f(freqs_p,Te,R)*exp(-1j*Ns*phi)

    bessels = _besselJ(Ns,z)

    return 0.5*(exp(-1j*z*_np.sin(phi))*_np.dot(bessels,Sds_m) + exp(+1j*z*_np.sin(phi))*_np.dot(bessels,Sds_p))
Sphi_of_f = _MemMut(_np.vectorize(_Sphi_of_f))

#@vectorize([float64(float64,float64,float64,float64,float64,float64,int32),float64(float64,float64,float64,float64,float64,float64,int64)])     
def Spa_adiabatique(omega, phi, nu, nuac, Te, R,phi_axis=-1) :
    """
    Voir note Udson Mendes
    phi_axis est l'axes des oscillations nuac
    """
    dphi  = phi[...,1] - phi[...,0]
    return (1./(2*_np.pi)) * Sdc_of_f( omega, nu + nuac , Te , R ).sum(axis=phi_axis) * dphi

#@vectorize([float64(float64,float64,float64,float64,float64,float64,int32),float64(float64,float64,float64,float64,float64,float64,int64)])     
def Spa_adia_harm_drive(omega,nu,nuac_amplitude,Te,R,n_phi=1000):
    """
    Adiabatique SII integrated over an harmonic drive
    """
    phi  = _np.linspace(0,2*_np.pi,n_phi)[:-1] # sinon le dernier point est en double
    Nuac = nuac_amplitude[:,None]*_np.sin(phi)[None,:]
    return Spa_adiabatique(omega[...,None],phi,nu[...,None],Nuac,Te[...,None],R[...,None],phi_axis=-1)

def _betap(p,f,nu,Te,nuac,Omega,R,nBessel=21):
    omega = 2.0*_np.pi*f
    if p==0:
        return Spa(omega,nu,Te,nuac,Omega,R,nBessel=nBessel)
    if not nuac*Omega:
        return 0
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
    freqs_m = -omega-Ns*Omega+nu
    freqs_p = +omega+Ns*Omega+nu

    Sds_m = Seq(freqs_m,Te,R)*(-1.)**p
    Sds_p = Seq(freqs_p,Te,R)

    bessels = _besselJ(Ns,z)*_besselJ(Ns+p,z)
    return _np.dot(bessels, (Sds_m+Sds_p)/2.)
betap = _MemMut(_np.vectorize(_betap))

def Xp(p,omega,nu,Te,nuac,Omega,R,nBessel=21):
    return betap(p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)+betap(-p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)

def Yp(p,omega,nu,Te,nuac,Omega,R,nBessel=21):
    return betap(p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)-betap(-p,omega,nu,Te,nuac,Omega,R,nBessel=nBessel)

def Sqz(p,omega,nu,Te,nuac,Omega, R,nBessel=21):
    return _np.sign(p)*betap(abs(p),omega,nu,Te,nuac,Omega,R,nBessel=21)+Spa(omega,nu,Te,nuac,Omega,R,nBessel=21)
    
    
######################
# Photon experiments #
######################

def C2_of_f(omega,nu,nuac=0.,Omega=0.,Te=0.050,R=50.0,nBessel=21):
    if (nuac == 0.) or (Omega == 0.):
        return (R/(_C.hbar*omega))*Sdc_of_f(omega,nu,Te,R) * 0.5
    else :
        z = nuac/Omega
        nBessel -= nBessel%2-1    # Ensure it's odd
        Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
        Sum = 0.
        for n in Ns:
            Sum +=  0.5*_besselJ(n,z)*_besselJ(n,z)*( Seq_of_f(nu+omega+n*Omega,Te,R) + Seq_of_f(nu-omega-n*Omega,Te,R) ) 
        return (R/(_C.hbar*omega))*Sum.sum() * 0.5 # this is necessary 
        
def dC2ac_th(omega,nu,nuac=0.,Omega=0.,Te=0.050,R=50.0,nBessel=21):
    return C2_of_f(omega,nu,nuac,Omega,Te,R,nBessel) - C2_of_f(omega,nu,0.,Omega,Te,R,nBessel)

def n_th(omega,nu,nuac=0.,Omega=1.e9,Te=0.050,R=50.0,nBessel=21):
    return C2_of_f(omega,nu,nuac,Omega,Te,R,nBessel) - 0.5 

def dnac_th(omega,nu,nuac=0.0,Omega=1e9,Te=0.050,R=50.0,nBessel=21):
    return n_th(omega,nu,nuac,Omega,Te,R,nBessel) - n_th(omega,nu,0.,Omega,Te,R,nBessel)

def C4_th(omega,nu,nuac,Omega,Te=0.050,R=50.0,nBessel=21):
    if omega == Omega/2 : # dangerous
        z = nuac/Omega
        nBessel -= nBessel%2-1    # Ensure it's odd
        Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
        Sum = 0.
        for n in Ns:
            Sum +=  0.5*_besselJ(n,z)*_besselJ(n+1,z)*( Seq_of_f(nu+omega+n*Omega,Te,R) - Seq_of_f(omega-n*Omega-nu,Te,R) ) 
        return ((3./2.)*((R/(_C.hbar*omega))*Sum.sum())**2) * 0.25 # this is necessary 
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
    return dn2_th(omega,nu,nuac,Omega,Te,R,nBessel)/n_th(omega,nu,nuac,Omega,Te,R,nBessel)

def Fanoac_th(omega,nu,nuac,Omega,Te=0.050,R=50.0,nBessel=21):
    return dn2ac_th(omega,nu,nuac,Omega,Te,R,nBessel)/nac_th(omega,nu,nuac,Omega,Te,R,nBessel)