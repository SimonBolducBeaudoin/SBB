#!/bin/env/python
#! -*- coding: utf-8 -*-

import os

def Copy_script(path_original,path_copy):
    """
    Missing reference...
    example :
    path_original = C:\\Projets\\Time_Domain_Optic\\Time_domain\\Python\\5-Experiments\\PQsVsVDC\\PQsvsVDC.py
    path_copy = C:\\Projets\\Time_Domain_Optic\\Time_domain\\Python\\5-Experiments\\TEST\\PQsvsVDC.py
    """
    with open(path_original, 'r') as f:
        with open(path_copy, 'w') as out:
            for line in (f.readlines()):
                out.write(line)
                    
def save_all_scripts(paths,scripts) :
    for key,value in scripts.items():
        Copy_script( paths['scripts']+os.sep+value , paths['saves']+os.sep+key+'.py')