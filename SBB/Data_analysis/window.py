#!/bin/env/python
#! -*- coding: utf-8 -*-

from __future__ import division
from past.utils import old_div
import numpy as _np

def window_after(S,i=65,t_demi=10):
    s = S.copy()
    epsilon = old_div(_np.log(2),t_demi)
    x    = _np.arange( S.shape[-1]-i )
    damp = _np.exp((-1)*epsilon*x)
    s[...,i:] *= damp
    return s