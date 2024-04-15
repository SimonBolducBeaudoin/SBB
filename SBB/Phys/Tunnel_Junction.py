#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np
from scipy.special import jv as _besselJ
import scipy.constants as _C

from SBB.Brrrr.Memoize import MemoizeMutable as _MemMut
from SBB.Numpy_extra.numpy_extra import find_nearest_A_to_a

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
    return _besselJ(n,z)**2.

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


def SII_of_f(freq,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    """
    Inputs 
    freq : 1D-array [Hz]                    (last axis of the final array)
    beta : ND-array or scalar [Hz^{-1/2}]   (last axis must be freq)
    
    Returns
    SII  : auto-correlation [A**2/Hz]       
    """
    omega = 2*_np.pi*freq
    Omega = 2*_np.pi*F
    nu    = _C.e*Idc*R/_C.hbar
    nuac  = _C.e*Iac*R/_C.hbar
    
    if ( ( type(Iac) != _np.ndarray) and ( Iac==0.0 ) ) or ( ( type(F) != _np.ndarray) and ( F==0.0 ) ) :
        return Sdc_of_f(omega,nu,Te,R)
    else :
        z = nuac/Omega
        nBessel -= nBessel%2-1    # Ensure it's odd
        Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
        sum_shape =  ( len(Ns), ) + _np.broadcast(Idc,Iac,freq,F).shape
        Sum = _np.zeros( sum_shape  )
        for i,n in enumerate(Ns):
            Sum[i,...] +=  0.5*_besselJ(n,z)*_besselJ(n,z)*( Seq_of_f(nu+omega+n*Omega,Te,R) + Seq_of_f(nu-omega-n*Omega,Te,R) ) 
        return Sum.sum(axis=0)  
        
def _Spa_two_freq(omega,nu,nuac,Omega,Te,R,nBessel):
    """
    Voir : Forgues et al. Experimental violation of Bell-like Inequalities By Electronic Shot Noise
    
    f + f2 = F : f est la plus petite des deux fréquences tel que f + f2 = F 
    <I(f1)I(f2)> = Sum_{n in {-inf,inf}} 0.5*alpha_n [Seq(f_n+) - Seq(f_n-)]
    f_npm = f+ n F pm nudc
    alpha_n = J_n( nuac/F ) J_{n+1}( nuac/F )
    
    nudc = eVac / h
    nuac = eVac / h
    
    """
    
    z = nuac/Omega
    nBessel -= nBessel%2-1    # Ensure it's odd
    Ns = _np.arange(-nBessel//2+1,nBessel//2+1)
    Sum = 0.
    for n in Ns:
        alpha_n = _besselJ(n,z)*_besselJ(n+1,z)
        omega_p =  omega+n*Omega+nu
        omega_m =  omega+n*Omega-nu
        Sum +=  0.5*alpha_n*( Seq_of_f(omega_p,Te,R) - Seq_of_f(omega_m,Te,R) ) 
    return Sum.sum()
Spa_two_freq = _np.vectorize(_Spa_two_freq)

def SII_f1_f2(freq,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    """
    f + f2 = F : f est la plus petite des deux fréquences tel que f + f2 = F 
    """
    omega = 2*_np.pi*freq
    Omega = 2*_np.pi*F
    nu    = _C.e*Idc*R/_C.hbar
    nuac  = _C.e*Iac*R/_C.hbar
    return Spa_two_freq(omega,nu,nuac,Omega,Te,R,nBessel)
    
def C4_f1_f2(freq,beta,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    _,F_idx = find_nearest_A_to_a(F,freq)
    df = freq[1]-freq[0]
    # normalised to 1/2 of int_+ df now will be 1.
    B = _np.absolute(beta)
    # Creating B(F-f)
    B2 = _np.zeros(B.shape)
    B2[...,:F_idx+1:] = B[...,F_idx::-1] # From F (or 12GHz) to 0 backward  
    B2[...,F_idx+1:]  = B[...,:len(freq)-(F_idx+1):][...,::-1] # The part of the spectrum folded around 0.
    
    sii_f1f2 = SII_f1_f2(freq,Idc,Iac,F,Te,R,nBessel) # |<i(f)i(F-f)>| [A^2/Hz]
    sii_f1f2[...,F_over2_idx] = 0.5 * sii_f1f2[...,F_over2_idx]
    
    vec1 = ((R/_C.h))*sii_f1f2/_np.sqrt(_np.abs(freq*(F-freq))) # has dimension of n2 but is not n2
    # Managing division by 0.
    w = _np.where( (vec1 == _np.inf) | (vec1 == -_np.inf)  ) 
    vec1[w] = 0.0
    
    vec2 = _np.zeros(vec1.shape)
    vec2[...,:F_idx+1:] = vec1[...,F_idx::-1] # From F (or 12GHz) to 0 backward  
    vec2[...,F_idx+1:]  = vec1[...,:len(freq)-(F_idx+1):][...,::-1] # The part of the spectrum folded around 0.
    
    vec1 = B*B2*vec1
    vec2 = B*B2*vec2
    
    #return (3./4.)*( 0.5*_np.nansum( vec1 ,axis=-1)**2 + _np.nansum( vec1**2 ,axis=-1) + _np.nansum( vec1*vec2 ,axis=-1) ) * df**2
    return (3./4.)*( 0.5*_np.nansum( vec1 ,axis=-1)**2 ) * df**2
    
def C4_f1_f2_old(freq,beta,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    """
    On integre le correlateur <i(f1)i(f2)>
    sur le domaine f1 + f2 = F en pondérant pour les beta correctements.
    Il s'agit d'un domaine en X, il faut éviter de compter le centre du X deux fois (i.e. f = F/2).
    
    F has to be in freq, because I need to find its index and use it in the computation.
    
    freq  = [0.2,0.4,...,F,F+0.2,...] 
    
    """
    _,F_idx = find_nearest_A_to_a(F,freq)
    _,F_over2_idx = find_nearest_A_to_a(freq[F_idx]/2.,freq)
    # normalised to 1/2 of int_+ df now will be 1.
    B = _np.absolute(beta)**2
    # Creating B(F-f)
    B2 = _np.zeros(B.shape)
    B2[...,:F_idx+1:] = B[...,F_idx::-1] # From F (or 12GHz) to 0 backward  
    B2[...,F_idx+1:]  = B[...,:len(freq)-(F_idx+1):][...,::-1] # The part of the spectrum folded around 0.
    
    #BB = B*_np.roll( B[...,::-1],F_idx+1,axis=-1 )
    c4_f1f2 = SII_f1_f2(freq,Idc,Iac,F,Te,R,nBessel)**2 # |<i(f)i(F-f)>|^2 [A^4/Hz^2]
    
    # The f = F/2 term is counted in double in the calculation for SII_f1_f2, 
    # and since we squared it we need to divide this frequency by 4 
    # to get the right amplitude
                               # 0.25 ??? ou 0.0
    c4_f1f2[...,F_over2_idx] = 0.25 * c4_f1f2[...,F_over2_idx]
    
    n2_tilde_f1f2 = ((R/_C.h)**2)*c4_f1f2/_np.abs(freq*(F-freq)) # has dimension of n2 but is not n2
    
    df = freq[1]-freq[0]
    # Managing division by 0.
    w = _np.where( (n2_tilde_f1f2 == _np.inf) | (n2_tilde_f1f2 == -_np.inf)  ) 
    n2_tilde_f1f2[w] = 0.0
    
    # It is not clear to me where where it comes from but there is a remaning 
    # multiplying constant between this theory and the data. 
    # It is probably due to some folding spectrum issue.
    # Nevertheless I algomerated this error in K which teh details
    # can be figured out latter.
    
    K = 8
 
    return K*(3./2.)*_np.nansum( B*B2* n2_tilde_f1f2 ,axis=-1)* df**2
    
        
def n_beta(freq,beta,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    """
    Ergo. Hyp. ==> <n> = < \bar{n(t)} > ==>
    <n> = int_+ |beta(f)|^2 ( SII(f)/Zhf - 1/2 ) df
    
    Inputs 
    freq : 1D-array [Hz] (last axis of the final array)
    beta : ND-array or scalar [Hz^{-1/2}] (last axis must be freq)
    SII  : auto-correlation [A**2/Hz] (last axis must be freq)
    
    Returns
    n : array with shape = beta.shape[:-1]
    """
    
    SII = SII_of_f(freq,Idc,Iac,F,Te,R,nBessel)
    
    df = freq[1]-freq[0]
    
    # Sprectum folding factor, 0.5, included (i.e. SII theory are usually integrated over the full spectrum).
    n_of_f = 0.5*R*SII/(_C.h * freq ) - 1./2. 
    # Managing division by 0.
    w = _np.where( (n_of_f == _np.inf) | (n_of_f == -_np.inf)  ) 
    n_of_f[w] = 0.0
 
    B = 2*_np.absolute(beta)**2 # normalise to 1/2 of int_+ df
    return _np.nansum( B * n_of_f ,axis=-1)* df
    
def dn2_beta(freq,beta,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21):
    """
    
    """
    n = n_beta(freq,beta,Idc,Iac=0.,F=0.,Te=0.050,R=50.0,nBessel=21)
    if ( ( type(Iac) != _np.ndarray) and ( Iac==0.0 ) ) or ( ( type(F) != _np.ndarray) and ( F==0.0 ) ) :
        return n*(n+1.0)
    else :
        # Sprectum folding factor, 0.5, included (i.e. SII theory are usually integrated over the full spectrum).
        c4 = 0.25*C4_f1_f2(freq,Idc,Iac,F,Te,R,nBessel)
    
    

    
