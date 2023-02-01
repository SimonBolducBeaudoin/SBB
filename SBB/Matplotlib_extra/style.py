#!/bin/env/python
#! -*- coding: utf-8 -*-

color_list      = ['b', 'g', 'r', 'c','m','y','k']
linestyle_list  = ['-','--','-.',':']
marker_list     = ['*','+','x','o','.','D','s',',','v','^','<','>','1','2','3','4'] 


def gen_cycler(**options):
    return cycler(color=color_list[:4],linestyle=linestyle_list[:4],marker=marker_list[:4])
    
if __name__ == "__main__" :
    import matplotlib
    ratio =1.62
    width = 14*1.5
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]

    font = {'family' : 'None',
            'weight' : 'None',
            'size'   : 15}

    matplotlib.rc('font', size=20)

    rcParams['lines.linewidth'] = 3
    rcParams['lines.markersize'] = 15
    rcParams['lines.markeredgewidth'] = 2.5