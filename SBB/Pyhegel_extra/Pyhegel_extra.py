#!/bin/env/python
#! -*- coding: utf-8 -*-

from __future__ import division
from past.utils import old_div
__all__ = ["analyseur_de_spectre_scalaire","parralel_get"]
from SBB.Brrrr.Threading_extra import ThreadWithReturnValue as _Thread_rv

def analyseur_de_spectre_scalaire(data,dt,l_chunk=1024):
    """
    Converts raw signal to a scalar spectrum of lenght l_chunk
    
    Can be used on the Guzik's output
    """
    f = rfftfreq(l_chunk,dt)
    F = zeros(len(f))
    N_chunk = len(data)//l_chunk
    for i in range(N_chunk):
        F += abs(rfft(data[i*l_chunk:(i+1)*l_chunk]))
    return f,old_div(F,N_chunk)
    
