#!/bin/env/python
#! -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt

color_list      = ['b', 'g', 'r', 'c','m','y','k']
linestyle_list  = ['-','--','-.',':']
marker_list     = ['*','+','x','o','.','D','s',',','v','^','<','>','1','2','3','4'] 

def gen_cycler(**options):
    return cycler(color=color_list[:4],linestyle=linestyle_list[:4],marker=marker_list[:4])

def trigger_style_0():
    ratio =1.62
    width = 14*1.5
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    matplotlib.rcParams['font.size'] = 20 
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['lines.markersize'] = 15
    matplotlib.rcParams['lines.markeredgewidth'] = 2.5
    
def trigger_style_1():
    ratio =1.62
    width = 14*0.75
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    matplotlib.rcParams['font.size'] = 20 
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['lines.markersize'] = 15
    matplotlib.rcParams['lines.markeredgewidth'] = 2.5

def trigger_style_2():
    ratio =1.62
    width = 14*1.5
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    matplotlib.rcParams['font.size'] = 30 
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['lines.markersize'] = 15
    matplotlib.rcParams['lines.markeredgewidth'] = 2.5

def set_prl_rc_params(**kwargs):
    """
    Set Matplotlib rc params for compatibility with PRL's format conventions.
    """
    params = {
        'font.size' : 8,
        'font.family': 'serif',
        'font.serif': ['CMU Serif',], # The font must be installed
        #'text.usetex': True, # Better latex integration, only works if tex is installed
        'axes.labelsize': 6,
        'axes.titlesize': 10,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 5,
        #'legend.title_fontsize': 6, # No worky in python2
        'figure.titlesize': 12,
        'figure.dpi' : 600 ,
        'lines.linewidth': 1,
        'lines.markersize': 1.5,
        'lines.markeredgewidth' : 1.5 ,
        'patch.linewidth': 1,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'legend.frameon': False,
        'legend.loc': 'best',
        'axes.linewidth': 0.5,
        'xtick.top': True,
        'ytick.right': True,
        'axes.grid': False,
        'grid.linestyle': '--',
        'grid.alpha': 0.7,
    }
    params.update(**kwargs)

    plt.rcParams.update(params)

def trigger_PRL(w,h,**rc_kwargs):
    D = {'figure.figsize' : [w,h]}
    D.update(**rc_kwargs)
    set_prl_rc_params(**D)

def trigger_PRL_single_col(ratio =1.62,**rc_kwargs):
    width = 3+3./8
    trigger_PRL(width,width/ratio,**rc_kwargs)
    
def trigger_PRL_1_5_col(ratio =1.62,**rc_kwargs):
    width = (3+3./8)*1.5
    trigger_PRL(width,width/ratio,**rc_kwargs)
    
def trigger_PRL_2_col(ratio =1.62,**rc_kwargs):
    width = (3+3./8)*2.0
    trigger_PRL(width,width/ratio,**rc_kwargs)
 
def trigger_default_style():
    trigger_style_0()
    
def trigger_small_style():
    trigger_style_1()
    
def trigger_beamer_style():
    trigger_style_2()

    
if __name__ == "__main__" :
    trigger_my_default_style()