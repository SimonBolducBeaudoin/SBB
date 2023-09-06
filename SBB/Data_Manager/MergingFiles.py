#!/bin/env/python
#! -*- coding: utf-8 -*-

"""Regroups function used to combine and clean up data produced by experiments"""

import numpy as _np

def get_untreated_files(folder,pattern= 'exp_??????-??????.npz',ReDo=False,AddIgnore=True,LookForLonelyIgnore=True): 
    """ 
    Used to get a list of files that have not been treated yet (aka not flagged with ignore yet). 
         
    Basic use (ReDo=False,AddIgnore=False,LookForLonelyIgnore=False) 
        Looks inside folder for files matching "pattern",  
        also looks for files matching "pattern"+.ignore  
        and return a list of pattern files excluding those that have a corresponding .ignore sister 
        Exemple :  
        If folder content is : 
             'exp_000000-000000.npz' 
             'exp_000000-000001.npz' 
             'exp_000000-000001.npz.ignore' 
        then it returns 
             ['exp_000000-000000.npz'] 
              
    Options : 
        ReDo (False by default) 
            Doesn't filter for .ignore file and hence returns all files matching the pattern 
        AddIgnore (True by default) 
            Adds .ignore to files name (modifies files names) 
        LookForLonelyIgnore 
            Looks if some "pattern.ignore" files dont have a matching "parttern" file and add them to the list 
    """ 
    pattern_list = glob.glob(folder+os.sep+pattern) 
    ignore_list  = glob.glob(folder+os.sep+pattern+'.ignore') 
     
    # New_list with only experiments that have not already been processed     
    not_done_list = [ exp for exp in pattern_list if (exp + '.ignore' not in ignore_list) ] 
 
    # Adds .ignore to the latter 
    if AddIgnore : 
        for exp in not_done_list : 
            os.rename(exp,exp+'.ignore') 
 
    # Get the "parttern.ignore" files that don't have a corresponding "pattern" file 
    lonely_ignore_list = [exp.replace('.ignore', '') for exp in ignore_list if ( (exp.replace('.ignore', '') not in pattern_list) and (LookForLonelyIgnore) ) ] 
     
    if ReDo : 
        return pattern_list + lonely_ignore_list 
    else : 
        return not_done_list + lonely_ignore_list     
    
def Append_NPZ(file,compress=True,**data):
    """
    Notes if the archive is compressed apperenttly the only way out is to decrompress and compress.
    for a non compressed zip file the code could be modified to simply append.
    """
    with dict(_np.load(file,allow_pickle=True)) as f:
        f.update(data)
        if compress :
            _np.savez_compressed(file, **f)
        else :
            _np.savez(file, **f)
    