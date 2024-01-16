#!/bin/env/python
#! -*- coding: utf-8 -*-

import matplotlib

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

def trigger_PRL(w,h):
    matplotlib.rcParams['figure.figsize'] = [w,h]
    matplotlib.rcParams['figure.dpi'] = 600 
    matplotlib.rcParams['font.size'] = 15
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['lines.markersize'] = 15
    matplotlib.rcParams['lines.markeredgewidth'] = 2.5
  
def trigger_PRL_single_col():
    ratio =1.62
    width = 3+3./8
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    trigger_PRL(width,width/ratio)
    
def trigger_PRL_1_5_col():
    ratio =1.62
    width = (3+3./8)*1.5
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    trigger_PRL(width,width/ratio)
    
def trigger_PRL_2_col():
    ratio =1.62
    width = (3+3./8)*2
    matplotlib.rcParams['figure.figsize'] = [width,width/ratio]
    trigger_PRL(width,width/ratio)
 
def trigger_default_style():
    trigger_style_0()
    
def trigger_small_style():
    trigger_style_1()
    
def trigger_beamer_style():
    trigger_style_2()

    
if __name__ == "__main__" :
    trigger_my_default_style()