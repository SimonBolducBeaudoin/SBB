#!/bin/env/python
#! -*- coding: utf-8 -*-

import numpy
from SBB.Pyhegel_extra.Experiment import logger,Experiment

def combine_arrays(S_list,V_list):
    """
    To be deleted 
    """     
    # Check that V_list like a list of list of 1D arrays
    if isinstance( V_list[0][0], (list, tuple, ndarray) ) :
        pass # List of list ok
    elif isinstance( V_list[0], (list, tuple, ndarray) ) :
        # Single list of conditions
        V_list = [ [v] for v in V_list ]
    else :
        raise Exception("V_list has to be like a list of list of 1D np.arrays")
        
    # Check that Vs all have the same length
    if all( [ len(V_list[0]) != len(V) for V in V_list ] ) :
        raise Exception("len (V_0) != len(V_n)")
    
    # Get the number of dimensions of the experimental conditions
    n_dim_cdn = len(V_list[0])
    
    # A list of unique conditions for each experimental variables
    V =[ sort(unique(concatenate([v[n] for v in V_list]))) for n in range(n_dim_cdn) ]
    
    # A list of list of lenght for each experimental variables of each experiments
    l_list = [[len(v) for v in V_rep] for V_rep in V_list ]
    
    # A List of reshaped measurment with all experimental variables flatten in a single dimension
    S_list = [ Sn.reshape( (Sn.shape[0],) + (prod(ln),) + Sn.shape[1+n_dim_cdn:] ) for Sn,ln in zip(S_list,l_list) ]
    
    S = []
    # Iterate over the combined experimental conditions
    for i, (v_idx, v) in enumerate( super_enumerate(*V) ):
        # A list og repetitions of v 
        S_reps = reperition_of_v(S_list,V_list,v)
        # Concatenate into an np.array along the repetition axis
        S.append(concatenate( S_reps, axis=0))
    
    S = pad_arrays( *(S), axis =0  ) # Pad only the repetition axis
    # Concatenate the arrays in the list into a single array
    S = concatenate(S, axis=1)
    
    # Shape of the unique experimental variables
    l = tuple( len(v) for v in V )
    # Restoring shapes
    S.shape = (S.shape[0],) + l + S.shape[2:]

    # Cleaning up/Removing index along the first dimension for which S is all 0 
    S = remove_nan_subarrays(S)

    # Return the combined arrays
    return V, S


class logger_acq_and_compute(logger):
    """
        A logger with default Aquisition and Computing events
    """
    def __init__(self,time_estimates,*arg,**log_dict):
        if log_dict :
            super(logger_aqc_and_compute,self).__init__(time_estimates,log_dict)
        else :
            n_measures  = arg[0]
            l_Vdc       = arg[1]
            l_data      = arg[2]
            _aqc_and_compute_log       = \
            {
            'loop_sizes'    : ( n_measures , l_Vdc ),
            'events'        : 
                            [
                                "Acquisition : {:04.2F} [s]", 
                                "Computing : {:04.2F} [s] "
                            ],
            'rate'          : ( l_data*1.0e-9 ,"Rate : {:04.2F} [GSa/s] " )
            }
            super(logger_aqc_and_compute,self).__init__(time_estimates,_aqc_and_compute_log)
            
logger_aqc_and_compute = logger_acq_and_compute

class Three_points_polarisation(object):
    """
       This class embeds the logic associated with 3 points measurment/polarisation
       (Voltmeter)      (Vdc)
       |                  |  
       |                 (R_pol) == 1MOhm : polarisation resistance
       |                  | 
       ----------|---------
                 |
                (R_s) : sample resistance
                 |
                (R_grnd) : ground/parasitic resistance
                 |
                ---
                 -
        Todos : 
            - Add imperfections in the junciton parasitic resistance capacitance ect.
            - Add photoexcitation behavior
        Bugs :
    """
    def __init__(self,R_s,R_pol):
        self.R_s    = R_s
        self.R_pol  = R_pol
        # self.R_grnd = R_grnd
    @staticmethod
    def compute_V_sample(Vdc,R_s,R_pol):
        return Vdc*R_s/(R_pol+R_s)
    @staticmethod
    def compute_V_yoko(V_s,R_s,R_pol):
        return V_s*(R_pol+R_s)/R_s
    @staticmethod
    def compute_I_sample(Vdc,R_pol):
        return Vdc/R_pol

class Conditions_logic(object):
    """
        THIS CLASS IS DEPRECATED AND IS LEFT HERE ONLY TO BE COMPATIBLE WITH OLD CODE.
        IT IS GOING TO BE AVENTUALLY REMOVED
        
        STATIC ONLY CLASS
        
        Methods with first argument beiing self must be called like
        Conditions_logic.method(object,*args,**kwargs)
        
        Manages the logic arround the condition tuples and the
        experimental conditions
        
        It adds some options attributes via the _set_options function
        and some Vdc attributes via the build_attributes function 
        
        Add the folowing options :
            Interlacing     : reference condition in between each conditions
            Vdc_antisym : anntisymetrize Vdc  
        Todos : 
            - Add the options of interlacing other than every other point
                ex : every 2nd point ... ref cnd cnd        ref cnd cnd
                or   every 3rd point ... ref cnd cnd cnd    ref cnd cnd cnd
        Bugs :
    """
    @staticmethod
    def _set_options(self,**options):
        self._conditions_options =   {'antisym':options.get('Vdc_antisym') }
        self._ref_options        =   {'interlacing': options.get('interlacing') , 'no_ref':options.get('no_ref')}
    @staticmethod
    def _build_attributes(self):
        # todos add options for the different scenarios of experiments
        # Vdc only experiment
        Vdc                     = self.Vdc
        conditions_options      = self._conditions_options
        ref_options             = self._ref_options
        self._Vdc_antisym       = Conditions_logic.add_antisym            (Vdc               ,**conditions_options )
        self._Vdc_exp           = Conditions_logic.add_ref_conditions     (self._Vdc_antisym ,**ref_options )
        self._conditions_exp    = self._Vdc_exp
        # Vdc and temperature experiment
        # ....
    #############
    # Utilities #
    #############
    @staticmethod
    def add_antisym(Vdc,**sym_options):
        return numpy.concatenate(([(-1.0)*Vdc[::-1],Vdc])) if sym_options.get('antisym') else Vdc
    @staticmethod
    def add_ref_conditions(Vdc,**ref_options):
        if    ref_options.get('no_ref'): 
            return Vdc
        elif  ref_options.get('interlacing'):
            return Conditions_logic.compute_interlacing(Vdc)
        else :
            return Conditions_logic.compute_default_ref(Vdc)
    @staticmethod
    def compute_interlacing(Vdc):
        Vdc_interlaced = numpy.zeros(2*len(Vdc))
        Vdc_interlaced[1::2] = Vdc
        return Vdc_interlaced
    @staticmethod
    def compute_default_ref(Vdc):
        return numpy.concatenate(([0],Vdc))
    #################
    # Loop behavior #
    #################
    @staticmethod
    def _core_loop_iterator(self):
        return Experiment._super_enumerate(self._conditions_exp)
    ######################
    # Analysis Utilities #
    ######################
    @staticmethod
    def get_conditions_slice(**ref_options):
        if ref_options.get('no_ref'):
            return slice(None)
        elif ref_options.get('interlacing'):
            return slice(1,None,2)
        else :
            return slice(1,None,None)
    @staticmethod
    def get_references_slice(**ref_options):
        if ref_options.get('no_ref'):
            raise ConditionsError('No references')
        elif ref_options.get('interlacing'):
            return slice(0,None,2)
        else :
            return [0,]
    @staticmethod
    def compute_cumulants_sample(cumulants,swapaxes=None,**ref_options):
        if swapaxes : 
            cumulants = cumulants.swapaxes(*swapaxes)
        slice_r             = Conditions_logic.get_references_slice(**ref_options)
        slice_c             = Conditions_logic.get_conditions_slice(**ref_options)
        ref                 = cumulants[...,slice_r]
        cdn                 = cumulants[...,slice_c]
        cumulants_sample    = cdn - ref
        if swapaxes : 
            cumulants_sample = cumulants_sample.swapaxes(*swapaxes)
        return cumulants_sample
    @staticmethod
    def compute_cumulants_sample_std(cumulants_std,swapaxes=None,**ref_options):
        if swapaxes : 
            cumulants_std = cumulants_std.swapaxes(*swapaxes)
        slice_r             = Conditions_logic.get_references_slice(**ref_options)
        slice_c             = Conditions_logic.get_conditions_slice(**ref_options)
        ref                 = cumulants_std[...,slice_r]
        cdn                 = cumulants_std[...,slice_c]
        cumulants_sample_std    = cdn + ref
        if swapaxes : 
            cumulants_sample_std = cumulants_sample_std.swapaxes(*swapaxes)
        return cumulants_sample_std

def _gen_dict_helper(d):
    out = dict()
    for k,i in list(d.items()):
        if i is not None:
            out.update({k:i})
    return out

def gen_quads_info(l_kernel,kernel_conf=None,alpha=None,filters_info=None):
    l_kernel    = int(l_kernel)
    l_hc        = l_kernel/2 + 1 
    quads_info  = {'l_kernel':l_kernel,'l_hc':l_hc,'kernel_conf':kernel_conf,'alpha':alpha,'filters_info':filters_info}
    return _gen_dict_helper(quads_info)
 
def gen_aqc_info(l_data,dt,gz_gain_dB=None,mv_per_bin=None):
    return _gen_dict_helper({'l_data':l_data,'dt':dt,'gz_gain_dB':gz_gain_dB,'mv_per_bin':mv_per_bin})
    
def gen_compute_info(n_threads,l_fft=None):
    return _gen_dict_helper({'n_threads':n_threads,'l_fft':l_fft})
    
def gen_circuit_info(R_jct,R_1M,R_tl=None,g=None,V_th=None):
    return _gen_dict_helper({'R_jct':R_jct,'R_1M':R_1M,'R_tl':R_tl,'g':g,'V_th':V_th})
 
def gen_hist_info(nb_of_bin,max):
    return _gen_dict_helper({'nb_of_bin':nb_of_bin,'max':max})