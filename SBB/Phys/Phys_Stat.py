#!/bin/env/python
#! -*- coding: utf-8 -*-

import scipy.constants as _C
import numpy as _np

def BoseEinstein(f,T):
    return 1.0/(_np.exp(_C.h*f/(_C.k*T))-1.0)