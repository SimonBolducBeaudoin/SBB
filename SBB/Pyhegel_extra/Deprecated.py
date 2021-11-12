#!/bin/env/python
#! -*- coding: utf-8 -*-

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
    __version__     = { 'Three_points_polarisation'  : 0.1 }
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
    __version__     = { 'Conditions_logic'  : 0.6 }
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