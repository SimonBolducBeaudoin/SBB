#!/bin/env/python
#! -*- coding: utf-8 -*-

import  sys, os, platform
import time

def gen_exp_default_paths(exp_root="5-Experiments",exp_dir="Unkown",module="Local_Module"):   
    paths = \
    {
        '_experiments_root'     : exp_root,                         # Relative to python_scripts_root
        '_current_experiment'   : exp_dir,                          # Relative to experiments_root
        '_lib'                  : module,                           # Relative to current_experiment
        '_test_dir'             : "TEST",                           # Relative experiments_root
        '_saves'                : time.strftime('%y%m%d-%H%M'),     # Relative to pwd/current_experiment 
    }
    paths['python_scripts_root']= "C:\\Projets\\Time_Domain_Optic\\Python_2_7"
    paths['experiments_root']   = paths['python_scripts_root']  + os.sep + paths['_experiments_root']  
    paths['current_experiment'] = paths['experiments_root']     + os.sep + paths['_current_experiment']  
    paths['Local_Module']       = paths['current_experiment']   + os.sep + paths['_lib'] 
    paths['scripts']            = paths['experiments_root']     + os.sep + paths['_current_experiment'] 
    paths['test']               = paths['experiments_root']     + os.sep + paths['_test_dir']
    paths['saves'] = paths['test']
    return paths

def gen_exp_default_scripts():
    """
        Les scripts par default d'une experience
    """     
    scripts = \
    {
        'Aquisition'            : 'Aquisition.py',
    }
    return scripts
    
def set_exp_environment(paths,scripts=None,test=False):
    """
    See : gen_exp_default_paths
    """        
    if not(test) :
        paths['saves'] = paths['current_experiment']+ os.sep + paths['_saves']    
    else : 
        paths['saves'] = paths['test']
        
    if scripts is None :
        scripts = gen_exp_default_scripts()
    try :
        from pyHegel.commands import make_dir
        make_dir(paths['saves']) # Tell to pyhegel to make the saves directory and to save data there.
    except :
        pass
    if paths['Local_Module'] not in sys.path :
        sys.path.append(paths['Local_Module'])
    return scripts,paths  

def add_Cygwin_mingw_to_path(mingw_path = 'C:\\cygwin64\\usr\\x86_64-w64-mingw32\\sys-root\\mingw\\bin' ):
    if os.name == 'nt':
        if platform.python_version_tuple()[0]== '2' : #python 2
            if ( mingw_path not in os.environ['PATH'] ) :
                os.environ['PATH'] = mingw_path+os.path.pathsep+os.environ['PATH']
        else : # python 3      
            os.add_dll_directory("C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin")
        