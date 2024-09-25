#!/bin/env/python
#! -*- coding: utf-8 -*-

from __future__ import division
from past.utils import old_div
import scipy.constants as _C
import numpy as _np

def BoseEinstein(f,T):
    return 1.0/(_np.exp(old_div(_C.h*f,(_C.k*T)))-1.0)