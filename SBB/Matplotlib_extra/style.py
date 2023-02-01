#!/bin/env/python
#! -*- coding: utf-8 -*-

color_list      = ['b', 'g', 'r', 'c','m','y','k']
linestyle_list  = ['-','--','-.',':']
marker_list     = ['*','+','x','o','.','D','s',',','v','^','<','>','1','2','3','4'] 

def gen_cycler(**options):
    return cycler(color=color_list[:4],linestyle=linestyle_list[:4],marker=marker_list[:4])