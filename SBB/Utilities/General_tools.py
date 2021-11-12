#!/bin/env/python
#! -*- coding: utf-8 -*-

from    scipy.optimize import leastsq
import  scipy.odr as odr  
import  numpy

from SBB.Time_quadratures.time_quadratures           import TimeQuad_uint64_t

import scipy.constants as C
_e      = C.e # Coulomb
_h      = C.h # Joule second
_hbar   = C.hbar
_kb     = C.Boltzmann # J/K (joules per kelvin)
_c      = C.c # m/s

def BoseEinstein(f,T):
    return 1.0/(numpy.exp(_h*f/(_kb*T))-1.0)
    
def gamma_to_Z(gamma,Z0=50.0):
    """
        Converts reflexion coefficient to impedance
        gamma = ZL - Z0 /  ZL + Z0*
    """
    return Z0*(gamma+1.0)/(1.0-gamma)

####################
# User inputs Utilities #
####################

def yes_or_no(question,default_ans='y'):
    while "the answer is invalid":
        reply = str(raw_input(question+' (y/n):').encode('utf-8')).lower().strip() or default_ans
        if reply[:1] == 'y':
            return True
        elif reply[:1] == 'n':
            return False
        else :
            return True if default_ans=='y' else False

####################
# math Utilities #
####################

def dB_to_V2_over_V1(dB,R2,R1=50.0):
    """
    V2/V1  = 10**(dBm/20) sqrt(R2/R1)
    """
    return numpy.sqrt(R2/R1)*10.0**(dB/20.0)
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
    return numpy.sqrt(2.0*R2*0.001)*10.0**(dBm/20.0)

def lin_to_dB(x):
    return 10*numpy.log10(x)

def dB_to_lin(x):
    return 10.0**(x/10.0)
    
def fourier_transform(F,dt):
    """
    This gives and approximation of F(f) or fourier transform of F(t) if you prefer
    and is not directly the DFT of F(t)    
    """
    return dt*numpy.fft.rfft(F)
def ifourier_transform(F,dt,n):
    """
    See Also
    -------
        fourier_transform
    """
    irfft = numpy.fft.irfft
    shift = numpy.fft.fftshift
    return (1.0/dt)*shift(irfft(F,n=n))

def moment(hx,hs,exp,n_total,no_clip=True):
    tmp = ((hx**exp)*hs)
    tmp = tmp[...,1:-1] if no_clip else tmp
    return (tmp.sum(axis=-1))/n_total
   
def centered_moment(hx,hs,exp,n_total,no_clip=True):
    mu = moment(hx,hs,1,n_total,no_clip)
    tmp = (((hx[None,None,:]-mu[...,None])**exp)*hs)
    tmp = tmp[...,1:-1] if no_clip else tmp
    return (tmp.sum(axis=-1))/n_total    

def polyfit_above_th(X,Y,Xth,deg=1,Y_even=None,swap=None):
    """
    Polynomial fit of Y(X) for X>Xth.
    If Y_even = True ==> 
        Fit Y(X) : X<-Xth and X>Xth
    If swap = True 
        swap axes on Y before anything else
    Parameters
    ----------
    X     : array like, shape (M) or (N,M)
    Y     : array like, shape (N,M) or (M,N) then use swap = True
    Xth  : float or array_like, shape or (N)
    deg   : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(N,deg+1) or (2,N,deg+1) [Y_even==True]
    
    See Also
    --------
    numpy.polyfit
    """
    if swap is not None : Y = numpy.swapaxes(Y,0,1)
    N           = Y.shape[0]
    M           = Y.shape[1]
    Xth         = numpy.array([Xth]) if type(Xth) is float else numpy.array(Xth)
    xth         = numpy.zeros((N,))
    xth[...]    = Xth
    x           = numpy.zeros((N,M))
    x[...]      = X
    if Y_even :
        X_pos       = x>=xth[:,None]
        X_neg       = x<=-xth[:,None]
        P  = numpy.zeros((2,N,deg+1))
        for j,(xx,y,x_pos,x_neg) in enumerate(zip(x,Y,X_pos,X_neg)):
            P[0,j]   = numpy.polyfit(xx[x_pos],y[x_pos],deg)
            P[1,j]   = numpy.polyfit(xx[x_neg],y[x_neg],deg)
    else :
        X_pos       = x>=xth[:,None]
        P  = numpy.zeros((N,deg+1))
        for j,(xx,y,x_pos) in enumerate(zip(x,Y,X_pos)):
            P[j]   = numpy.polyfit(xx[x_pos],y[x_pos],deg)
    return P

class fit:
    
    '''
    Uses least square regression to fit a given model.
    
    Parameters
    ----------
    x : np array or array like
        abscisse
    y : np array or array like
        ordonnée/data
    p0 : np array or array like
        The starting estimate for the minimization
        As to behave like a 
    yerr : np array or array like
    weigths : np array or array like
        Weigth attributed for each point. Higher weight means that this data points is more important to fit 
        properly while a weight near 0 means that this data point is uselees
    
    
    Tutorial
    ----------
        How to fit multidimensionnal arrays ?
            
    
    Todos
    ------
        - Update __doc__
        - Make tutorial
        - Trouver une manière élégante de faire des paquets de fits
        - Trouver une manière élégante de faire des fits 2D
        - I'm not sur xerr and yerr where implemented correctly in odr
    Bugs
    ----
    
    Last update notes
    -----------------
        Removed unecessary self.xFit
    
    Old documentation
    ----------------
    Modifier à partir de la version de Max
     
    La classe 'fit' crée un objet contenant deux méthodes de fit, soient
    fit.leastsq et fit.odr qui renvoient toutes deux le résultat dans
    fit.para avec comme erreurs fit.err. Si l'option verbose est activée
    (activée par défaut), le fit imprime à l'écran des valeurs importantes
    (dont la matrice de corrélation et chi^2).
    
    Elle s'utilise comme suit :
    def Fonction(x,P):
        return     La_fonction_de_la_variable_x_et_du_tableau_de_paramètres_p

    a = fit(ValeursDeX, ValeursDeY, ParamètresInitiaux, Fonction, xerr=ErreursEnX, yerr=ErreursEnY)
    a.leastsq() OU a.odr()

    Aussi, appeler l'objet de 'fit' correspond à appeler la fonction avec
    les paramètres stockés dans fit.para (paramètres initiaux au paramètres
    de fit)
    a(x) est absolument équivalent à Fonction(x,a.para)

    Aussi, on peut aller chercher directement les paramètres du fit en
    considérant l'objet comme un tableau:
    a[i] est absolument équivalen.let à a.para[i]

    Les classes 'lsqfit' et 'odrfit' sont absolument identiques à la classe
    'fit' (elles héritent de toutes ses méthodes et variables), sauf qu'elle
    performent la régression au moment de l'initialisation. Ainsi :
    def Fonction(x,P):
        return La_fonction_de_la_variable_x_et_du_tableau_de_paramètres_p

    a = odrfit(ValeursDeX, ValeursDeY, ParamètresInitiaux, Fonction, xerr=ErreursEnX, yerr=ErreursEnY)
        
    '''
    __version__ = {'fit':0.3}
    def __init__(self,x,y,p0,f,yerr=None,weights=None,verbose=True,fullo=1):
        self.x      = numpy.array(x)
        self.y      = numpy.array(y).flatten()
        yerr        = yerr      if yerr     is not None else [1.]*self.x.shape[-1]
        yerr        = numpy.array(yerr).flatten()
        weights     = weights   if weights  is not None else [1.]*self.x.shape[-1]
        weights     = numpy.array(weights).flatten()
        self.yerr   = yerr/weights        
        self.para   = numpy.r_[p0].flatten() # Paramètres initiaux
        self.f       = f        # Fonction pour la régression
        self.fullo   = fullo     # 'Full output' de leastsq
        self.verbose = verbose    # Imprime des résultats importants à l'écran
    ###########
    # dunders #
    ###########
    def __call__(self,x=None,p=None):
        if x is None:
            x = self.x
        if p is None:
            p = self.para
        return self.f(x,p)
    def __getitem__(self,i):
        return self.para[i]
    def __len__(self):
        return len(self.para)
    def __str__(self):
        return '\n--- FIT ON FUNCTION {} ---\n\nFit parameters are {}\nFit errors are {}\nFit covariance\n{}\nFit correlation matrix\n{}\nReduced chi2 is {}\n\n'.format(self.f.__name__,self.para,self.err,self.cv,self.corrM,self.chi2r)
    ###########
    # Private #
    ###########
    def _residuals(self,p):
        return (self.y-self.f(self.x,p))/self.yerr
    def _computevalues(self):
        self.sdcv = numpy.sqrt(numpy.diag(self.cv))
        # Matrice de corrélation
        self.corrM = self.cv/self.sdcv/self.sdcv[:,None]
        self.chi2 = sum(((self.y-self.f(self.x,self.para))/self.yerr)**2.)
        # Chi^2 réduit
        self.chi2r = self.chi2/(len(self.y)-len(self.para))
    ########
    # User #
    ########
    def diff(self):
        """
            Return the difference between the fitted model and the data
        """
        return self.y - self.f(self.x,self.para)
    def leastsq(self):
        """
            Computes the least squarre 
        """
        self.lsq = leastsq(self._residuals,self.para,full_output=self.fullo)
        if self.lsq[1] is None:
            if self.verbose: print '\n --- FIT DID NOT CONVERGE ---\n'
            self.err = None
            self.chi2r = None
            return False
        else:
            # Paramètres :
            self.para = self.lsq[0]
            self.cv = self.lsq[1]
            # Nombre d'itérations :
            self.it = self.lsq[2]['nfev'] 
            self._computevalues()
            self.err = self.sdcv*numpy.sqrt(self.chi2r)
#            self.donnee = []
#            for i in range(len(self.para)):
#                self.donnee.append(d.donnee(self.para[i],self.err[i]))
            if self.verbose:
                print self
            return True
    def model(self,x,p):
        """
            Used to evaluated the model at with given p and x
        """
        return self.f(x,p)
class lsqfit(fit):
    def __init__(self,x,y,p0,f,fullo=1,yerr=None,weights=None,verbose=True):
        fit.__init__(self,x,y,p0,f,fullo=fullo,yerr=yerr,weights=weights,verbose=verbose)
        self.leastsq()    

####################
# Python tricks Utilities #
####################
   
def formated_tuple(frmt='{:0.2f}',T=tuple()):
    s = ''
    T = numpy.array([T]) if (type(T)==float or type(T) == numpy.float64) else T 
    for t in T :
        s += frmt.format(t)
        s += ', '
    return s
   
def build_array_of_objects(shape,constructor,*args,**kargs):
    A = numpy.r_[[constructor(*args,**kargs) for i in range(numpy.prod(shape))]]
    A.shape = shape
    return A

####################
# Numpy tricks #
####################

def get_index(Xs,x):
    """
    Returns the single index for wich Vdc = V
    Returns false if it doesn't exist
    """
    tmp = where(Vdc==V)[0]
    return False if tmp.size==0 else tmp[0]

def find_nearest_A_to_a(a,A):
    a = numpy.array([a]) if type(a)==float else numpy.array(a)
    A = numpy.array(A)
    a_shape = a.shape
    a       = a.flatten()
    X       = numpy.empty(a.shape)
    X_idxs  = numpy.empty(a.shape,dtype=int)
    for i,x in enumerate(a) :
        index       = numpy.abs(A-x).argmin()
        X_idxs[i]   = index 
        X[i]        = A[index]
    X.shape         = a_shape
    X_idxs.shape    = a_shape
    return X , X_idxs 

def symetrize(X):
    """
        Symetrize along the last axis
    """
    return numpy.concatenate((X[...,-1:0:-1],X),axis=-1)
        
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
    X_diff          = numpy.zeros(shape)
    X_diff[0,...]   = X[...,1:]
    X_diff[1,...]   = X[...,:-1]
    return X_diff
    
def cyclic_tansformation(X):
    """
        Cyclic translation of nd array
    """
    shape     = X.shape
    out       = numpy.zeros(X.size)
    flat_input = X.flatten()
    out[0:-1] = flat_input[1:]
    out[-1]   = flat_input[0]
    out.shape = shape
    return out

####################
# Matplotlib tricks Utilities #
####################    

color_list      = ['b', 'g', 'r', 'c','m','y','k']
linestyle_list  = ['-','--','-.',':']
marker_list     = ['*','+','x','o','.','D','s',',','v','^','<','>','1','2','3','4']  

def gen_cycler(**options):
    return cycler(color=color_list[:4],linestyle=linestyle_list[:4],marker=marker_list[:4])

"""
    See :
    Executing modules as scripts
    https://docs.python.org/3/tutorial/modules.html
"""
if __name__ == "__main__":
    import sys
    fib(int(sys.argv[1]))
  