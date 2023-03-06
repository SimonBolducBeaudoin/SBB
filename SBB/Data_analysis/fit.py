#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np

def polyfit_above_th(X,Y,Xth,deg=1):
    """
    Polynomial fit of Y(X) for X>Xth.

    Inputs
    ------
    X     : array like, shape (M) or (N,M)
    Y     : array like, shape (N,M) or (M,N) then use swap = True
    Xth  : float or array_like, shape or (N)
    deg   : int
        Degree of the fitting polynomial
    
    Returns
    -------
    P     : ndarray, shape(N,deg+1)
    
    See Also
    --------
    _np.polyfit
    """
    N           = Y.shape[0]
    M           = Y.shape[1]
    Xth         = _np.array([Xth]) if type(Xth) is float else _np.array(Xth)
    xth         = _np.zeros((N,))
    xth[...]    = Xth
    x           = _np.zeros((N,M))
    x[...]      = X
    X_pos       = x>=xth[:,None]
    P  = _np.zeros((N,deg+1))
    for j,(xx,y,x_pos) in enumerate(zip(x,Y,X_pos)):
        # Removing nans 
        not_nan = ~(_np.isnan(y)) 
        all_nan = all(not_nan==False)
        if all_nan :
            P[j] = nan
        else :
            P[j]   = _np.polyfit(xx[x_pos & not_nan],y[x_pos & not_nan],deg)
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
        self.lsq = leastsq(self._residuals,self.para,full_output=self.fullo)
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