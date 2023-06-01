#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np
from scipy.optimize import leastsq as _leastsq

from SBB.Numpy_extra.numpy_extra import sub_flatten

import pdb

from matplotlib.pyplot import subplots

def polyfit_multi_check(x,Y,Xth=None,deg=1,Y_idxs=None,P=None,ax=None):
    """
    Fast visuallization of a polyfit_multi 
    
    Inputs
    ------
    x     : array like, shape (    M) 
    Y     : array like, shape (...,M) 
    deg   : int
        Degree of the fitting polynomial
    Y_idxs : list of length = ... (aka len(Y_idxs) = Y.ndim-1 )
        The Y idexes to plot
        will be randomly selected by default
    P : The already fitted P 
        will be fitted if not given
    ax : The axis for the plot
    
    Returns
    -------
    ax
    
    See Also
    --------
    polyfit_multi
    
    By default it selects one random index of Y[...,M] (aka a random index in the ... part)
    """
    def get_rdm_idx(n):
        if n == 0 :
            return 0
        else :
            return _np.random.randint(0,n)
    
    if Y_idxs is None :
        # selects a random index to plot
        Y_idxs = [ get_rdm_idx(Y.shape[i]-1) for i in range(Y.ndim -1) ] 
    
    if P is None :
        if Xth is None :
            P = polyfit_multi(x,Y,deg=deg) 
        else :
            P = polyfit_above_th(x,Y,Xth,deg=deg)
    
    if ax is None :
        fig,ax = subplots(1,1)
    
    line, =  ax.plot(x,Y[tuple(Y_idxs)],marker='.',ls='None')
    fit = 0.0
    for i,p in enumerate(P[tuple(Y_idxs)][::-1]) :
        fit += p*x**i
    ax.plot(x,fit,marker='None',ls='--',label='{}'.format(Y_idxs),color=line.get_color())
    ax.legend(title='Y_idxs')
    
    return ax , P

def polyfit_multi(x,Y,deg=1):
    """
    Polynomial fit of Y(X) along the last axis 
    
    Bonus : it works when Y contains nans.

    Inputs
    ------
    X     : array like, shape (    M) 
    Y     : array like, shape (...,M) 
    deg   : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(...,deg+1)
    
    See Also
    --------
    _np.polyfit
    """
    if Y.ndim == 1 :
        Y = Y[None,:]
    
    P_shape = Y.shape[:-1]+(deg+1,)
    P       = _np.full( (_np.prod(P_shape[:-1]),) + P_shape[-1:] , _np.nan )
       
    for j,y in enumerate( sub_flatten(Y) ):  #coulb be made faster by usign sub_flat instead
        not_nan = ~(_np.isnan(y)) 
        all_nan = all( not_nan ==False)
        if all_nan :
            P[j] = _np.nan
        else :
            P[j]   = _np.polyfit(x[not_nan],y[not_nan],deg)
    P.shape = P_shape
    return P

def polyfit_above_th(x,Y,Xth,deg=1):
    """
    For retrocompatibility
    """
    return polyfit_multi_between(x,Y,Xth_low=Xth,Xth_high=None,deg=deg)

def polyfit_multi_between(x,Y,Xth_low=None,Xth_high=None,deg=1):
    """
    Polynomial fit of Y(X) for (X>=Xth_low)&((X=<Xth_high)).

    Inputs
    ------
    X        : array like,                 shape (    M) 
    Y        : array like,                 shape (...,M) 
    Xth_low  : float or array_like with     shape (...  )
    Xth_high : float or array_like with     shape (...  )
    deg      : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(...,deg+1)
    
    See Also
    --------
    _np.polyfit
    """
    if (Xth_low is None) and (Xth_high is None):
        return polyfit_multi(x,Y,deg=deg)
    
    def build_Xth(Xth,Y):
        xth         = _np.full( Y.shape[:-1],_np.nan )
        xth[...]    = Xth
        return xth
        
    def select_x(x,xth_l,xth_h):
        if _np.isnan(xth_l) :
            return x<=xth_h
        elif _np.isnan(xth_h) :
            return x>=xth_l
        else :
            return (x>=xth_l)&(x<=xth_h)
    
    if Y.ndim == 1 :
        Y = Y[None,:]   
    Xth_low  = build_Xth(Xth_low ,Y)  # Memory hungry ... 
    Xth_high = build_Xth(Xth_high,Y)
    
    P_shape = Y.shape[:-1]+(deg+1,)
    P       = _np.full( (_np.prod(P_shape[:-1]),) + P_shape[-1:] , _np.nan )
       
    for j,(y,xth_l,xth_h) in enumerate(zip(sub_flatten(Y),Xth_low.flat,Xth_high.flat)):  #coulb be made faster by usign sub_flat instead
        x_pos = select_x(x,xth_l,xth_h) # Select the abscisse idx repecting the threshold
        not_nan = ~(_np.isnan(y))       
        all_nan = _np.all( (x_pos & not_nan) ==False)
        if all_nan :
            P[j] = _np.nan
        else :
            P[j]   = _np.polyfit(x[x_pos & not_nan],y[x_pos & not_nan],deg)
    P.shape = P_shape
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
        self.x      = _np.array(x)
        self.y      = _np.array(y).flatten()
        yerr        = yerr      if yerr     is not None else [1.]*self.x.shape[-1]
        yerr        = _np.array(yerr).flatten()
        weights     = weights   if weights  is not None else [1.]*self.x.shape[-1]
        weights     = _np.array(weights).flatten()
        self.yerr   = yerr/weights        
        self.para   = _np.r_[p0].flatten() # Paramètres initiaux
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
        self.sdcv = _np.sqrt(_np.diag(self.cv))
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
        self.lsq = _leastsq(self._residuals,self.para,full_output=self.fullo)
        if self.lsq[1] is None:
            if self.verbose: print('\n --- FIT DID NOT CONVERGE ---\n')
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
            self.err = self.sdcv*_np.sqrt(self.chi2r)
#            self.donnee = []
#            for i in range(len(self.para)):
#                self.donnee.append(d.donnee(self.para[i],self.err[i]))
            if self.verbose:
                print(self)
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