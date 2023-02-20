#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np

def window_after(S,i=65,t_demi=10):
    s = S.copy()
    epsilon = _np.log(2)/t_demi
    x    = _np.arange( S.shape[-1]-i )
    damp = _np.exp((-1)*epsilon*x)
    s[...,i:] *= damp
    return s