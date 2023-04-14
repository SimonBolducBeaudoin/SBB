#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy as _np
from scipy import constants  as C

def BoseEinstein(f,T):
    return 1.0/(_np.exp(C.h*f/(C.k*T))-1.0)