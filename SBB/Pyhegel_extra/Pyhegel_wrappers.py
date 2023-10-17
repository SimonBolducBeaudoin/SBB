#!/bin/env/python
#! -*- coding: utf-8 -*-

"""
    Pyhegel_wrapper should only work in a pyHegel instance
"""

# Pucblic 
__all__ = ["PyhegelWrapperErrors","CalibrationTableError","Pyhegel_wrapper","Dummy","Lakeshore_wrapper","Dmm_wrapper","Yoko_wrapper","Guzik_wrapper","PSG_wrapper"]


#Private
import numpy as _np , time as _time
import matplotlib.pyplot as _plt

from SBB.Histograms.histograms import Histogram_uint64_t_int16_t as Hist_int16_t , Histogram_uint64_t_double as _Hist_double
from SBB.Numpy_extra.numpy_extra import sub_flatten as _sub_flatten, build_array_of_objects as _build_array_of_objects

try :
    # pyHegel will not load properly if Qt is not active
    from pyHegel import instruments
    from pyHegel.commands import set,get, wait
except ImportError :
    pass
    
class PyhegelWrapperErrors(Exception):
    """
        Base class for pyHegel tools errors
    """
    pass
    
class CalibrationTableError(PyhegelWrapperErrors):
    """
        Called when outside the calibrated table's range
    """
    pass

####################
# Pyhegel wrapeprs #
####################

class Pyhegel_wrapper(object):
    """ 
        default values for every mandatory methods
    """
    __version__ = {'Pyhegel_wrapper':0.2}
    def get_pyhegel_instance(self):
        pass
    def set_init_state(self,val):
        pass
    def set_close_state(self):
        pass
    def get(self):
        pass
    def set(self,val):
        pass
             
class Dummy(Pyhegel_wrapper):
    """
        A dummy device that implement only default behavior and prints what it does.
    """
    def get_pyhegel_instance(self):
        print("get_pyhegel_instance")
    def set_init_state(self,val):
        print("set_init_state :{}".format(val))
    def set_close_state(self):
        print("set_close_state")
    def get(self):
        print("get")
    def set(self,val):
        print("set :{}".format(val))

class Lakeshore_wrapper(Pyhegel_wrapper):
    """ 
        Lakeshore_wrapper is a logical device derived from lakeshore_370
        
        Warning : this is just what works for me. I did not code this denfensivly
        use at you're own risk.
        
        It is intended to automatically control the temparature using the lakeshore's internal 
        closed loop along with presets (static tables) that I've made myself.
        
        Usgin object.stabilize_temperature(t) will give back the control of to the command line
        only when the temperature is reached and has stabilized (see : time_tolerance and temperature_tolerance)
        
        If temperatures out of the intended range are used it defaults to 
        cooldown (i.e. self.heater_range is set to 0).
        
        Notes :
            - I changes getout_tunder() heating_power from 12.0e-6 to 6.0e-6 because 
              it was overshooting to much.
        Todos :
            - Option to show the power used during th
        Bugs :
            - Lorsque cooldown a ete appeler (fin du script) et que je relance le même script
            le lakeshore ne semble pas chauffé. J'ai du fermer le heater et le réouvrir en mettant un range pour que ça fonctionne comme d'hab.
    """
    __version__ = {'Lakeshore_wrapper':0.5}
    __version__.update(Pyhegel_wrapper.__version__)
    # These tables could be fine tuned
    # And are calibrated on  dilu 3
    _temperature_range       =   [ 0.00700005 , 0.015   , 0.020     , 0.050     , 0.100    , 0.200     , 0.400  , 0.700  , 1.1  ] # The upperbound of the interval
    _temperature_tolerance   =   [ 0.0      , 0.0005    , 0.00025   , 0.0005    , 0.00025  , 0.001     , 0.005  , 0.010  , 0.020] # The upperbound of temperature error after _stabilization_loop
    _time_tolerance          =   [ 120      , 600*3     , 222*3     , 105*3     , 250*3    , 420*3     , 600*3  , 600*3  , 600*3] # The amount of time the temperature error as to be below the temperature_tolerance to say the stabvilization is acheived
    _heater_range            =   [ 0.0      , 0.000316  , 0.000316  , 0.001     , 0.00316  , 0.01      , 0.01   , 0.01   , 0.01 ] # see instruments.lakeshore_370.heater_range for all possible values
    _proportional            =   [ 25       , 25        , 25        , 25        , 12.5     , 6.25      , 6.25   , 6.25   , 6.25 ]
    _integral                =   [ 120      , 120       , 120       , 120       , 120      , 120       , 120    , 120    , 120  ]
    _derivative              =   [ 120        , 0       , 0         , 0         , 0        , 0         , 0      , 0      , 0    ]
    _observed_time_cst       =   [ 1600     , 600       , 220       , 105       , 250      , 420       , 600    , 600    , None ] # [s] very approximate
    _base_flow               =   0.50                                                                                     # mmol/s Base flow on this fridge for this calibration table
    def __init__(self, visa_addr ='ASRL1',**options):
        self._set_options(**options)
        if self._debug :
            self._lakeshore = None
        else :
            self._lakeshore =  instruments.lakeshore_370(visa_addr)
            self.set_current_ch(6) # The 6th channel is the temperature of the sample
        if self._verbose :
            print("Current channel is : {} \n and has temperature = {} K".format( self._lakeshore.current_ch.get(), self.get_temperature()))
    #############
    # User interface #
    #############
    def get_pyhegel_instance(self):
        return self._lakeshore
    def _getout_tunder(self,verbose=None):
        """
            Runs the gets out t.under routine only if necessary
            Todos :
                - 
            Bugs :
                - Verbose not showing in command line
        """
        if self._debug :
            pass
        else :
            if (self._verbose or verbose) :
                print("get out of t.under")
            ch_sample       = 6         # ch=6 : mixing chamber
            T_target_out    = 0.015     # = 15[mK] : get the value from table 
            heating_power   = 6.0e-6    # [W] 100% of ther targeted range   
            t_th            = self._temperature_range[0]    # Thereshold to say we are out of t.under 
            wait_time       = 5.0       # Before checking if we're out of t.under
            if self.get_temperature() <= t_th :
                self._set_autoscan(False)
                self._set_autoscan_ch(ch_sample)    
                self._set_to_open_loop()
                self._set_heater_range(T_target_out)    
                self._set_heater(heating_power)
                while True :
                    t = self.get_temperature()
                    if t > t_th :
                        break
                    if (self._verbose or verbose) :
                        print("Waiting to get out of t.under")
                    wait(wait_time)
            if (self._verbose or verbose) :
                print("setting lakeshore to close loop")
            self._set_to_close_loop()
    def cooldown(self):
        if self._debug :
            pass
        else :
            self._set_heater_range(0)    # Need to turn the heater off first
            self._set_control_off()      # Then set the control mode to off (no temperature control)
    def set_close_state(self):
        if self._debug :
            pass
        else :
            self.cooldown()
    def set_current_ch(self,channel):
        if self._debug :
            pass
        else :
            set(self._lakeshore.current_ch,channel)
    def get_current_ch(self):
        if self._debug :
            return None
        else :
            return get(self._lakeshore.current_ch)
    def get_temperature(self):
        if self._debug :
            return 0.0
        else :
            return get(self._lakeshore.t)
    # def get(self):
        # return self.get_temperature()
    # def set(self,T):
        # return self.stabilize_temperature(T)
    def set_temperature(self,t):
        if self._debug :
            pass
        else :
            set(self._lakeshore.sp,t)
    def stabilize_temperature(self,T,verbose=False):
        if self._debug :
            pass
        else :
            if (self._verbose or verbose) :
                print("Stabilizing to {}[K]".format(T))
            try :
                if T > self._temperature_range[0] :
                    self._getout_tunder(verbose)
                self._set_heater_range(T)
                self._set_PID(T)
                self.set_temperature(T)
                self._stabilization_loop(T)
            except CalibrationTableError :
                self.cooldown()
                raise CalibrationTableError(' Target temperature was above table''s range. \n Heater was set to 0 and control_mode to off')
    #############
    # private #
    #############
    def _set_options(self,**options):
        self._debug     = options.get('debug')
        self._verbose   = options.get('verbose')
    def _set_autoscan(self,val_bool):
        set(self._lakeshore.scan,autoscan_en=val_bool)
    def _set_autoscan_ch(self,ch):
        set(self._lakeshore.scan,ch=ch)
    def _set_control_off(self):
        set(self._lakeshore.control_mode,'off') 
    def _set_to_open_loop(self): 
        set(self._lakeshore.control_mode,'open_loop')
    def _set_to_close_loop(self):
        if get(self._lakeshore.control_mode)!='pid':
            set(self._lakeshore.control_mode,'pid')
    def _get_table_index(self,T):
        t_table = self._temperature_range
        if not(  T<=t_table[-1] ) :
            raise CalibrationTableError   
        for i , t in enumerate(t_table):
            if T <= t :
                index = i  
                break 
        return index  
    def _stabilization_loop(self, T_target):
        t_table = self._temperature_range
        if T_target <= t_table[0] :
            T_target = t_table[0]
        t_tol       = self._temperature_tolerance[self._get_table_index(T_target)]
        time_tol    = self._time_tolerance[self._get_table_index(T_target)]
        T           = self.get_temperature()
        delta_T     = _np.abs(T - T_target)
        keep_going  = True
        converged   = False
        while keep_going :
            t_0 = _time.time()
            t_1 = t_0
            while delta_T<=t_tol :
                wait(5)
                T       = self.get_temperature() 
                delta_T = _np.abs(T - T_target)
                t_1 = _time.time()
                if not delta_T<=t_tol :
                    break 
                elif (t_1-t_0)>=time_tol :
                    converged = True
                    break 
            if converged == True : 
                keep_going = False
                break
            wait(10)
            T       = self.get_temperature() 
            delta_T = abs(T - T_target)       
    def _set_heater(self,power):
        """ Power in [W]"""
        set(self._lakeshore.manual_out_raw,power)
    def _set_heater_range(self, T):
        index = self._get_table_index(T)
        set(self._lakeshore.heater_range,self._heater_range[index])
    def _set_PID(self,T):
        index = self._get_table_index(T)
        set(self._lakeshore.pid, P = self._proportional[index] )
        set(self._lakeshore.pid, I = self._integral[index]     )
        set(self._lakeshore.pid, D = self._derivative[index]   )

class Dmm_wrapper(Pyhegel_wrapper):
    """
    """
    __version__ = {'Dmm_wrapper':0.1}
    __version__.update(Pyhegel_wrapper.__version__)
    name_addr = {  'dmm11':'GPIB0::22::INSTR', }
    def __init__(self,visa_addr=None,nickname='dmm11',**options):
        self._set_options(**options)
        if self._debug :
            self._dmm = None
            self._V_dummy = 0.0
        else :
            visa_addr = visa_addr if visa_addr else Dmm_wrapper.name_addr[nickname]
            self._dmm = instruments.agilent_multi_34410A(visa_addr)
            set(self._dmm.aperture, 1)
            set(self._dmm.zero,True)
    """
        Wrapper of existing behavior
    """
    def _set_options(self,**options):
        self._debug  = options.get('debug')
    def get(self):
        if self._debug :
            print "get Dmm : {:0.2f}[V]".format(self._V_dummy)
            return self._V_dummy
        else :
            return get(self._dmm)

    def set_aperture(self,val=None):
        if self._debug :
            print "set Dmm aperture : {}".format(val)
        else :
            set(self._dmm.aperture, val)
    def set_zero(self,booleen=True):
        if self._debug :
            print "set Dmm zero : {}".format(booleen)
        else :
            set(self._dmm.zero, booleen)
        
class Yoko_wrapper(Pyhegel_wrapper):
    """
        Todos :
            - Add debug option 
    """
    __version__ = {'Yoko_wrapper':0.4}
    __version__.update(Pyhegel_wrapper.__version__)
    name_addr = {  'yo10':'USB0::0x0B21::0x0039::91KB11655',
                   'yo11':'USB0::0x0B21::0x0039::91M504010'
                }
    def __init__(self,visa_addr=None,nickname='yo10',**options):
        self._set_options(**options)
        if self._debug :
            self._yoko = None
            self._V_dummy = 0.0
        else :
            visa_addr = visa_addr if visa_addr else Yoko_wrapper.name_addr[nickname]
            self._yoko = instruments.yokogawa_gs200(visa_addr)
        self.set(0)
    """
        Wrapper of existing behavior
    """
    def _set_options(self,**options):
        self._debug  = options.get('debug')
    def get(self):
        if self._debug :
            print("get Yoko : {:0.2f}[V]".format(self._V_dummy))
            return self._V_dummy
        else :
            return get(self._yoko)
    def set(self,V):
        if self._debug :
            self._V_dummy = V
            print("set Yoko : {:0.2f}[V]".format(V))
        else :
            set(self._yoko,V)
    def set_output(self,booleen):
        if self._debug :
            print("set Yoko outpout : {}".format(booleen))
        else :
            set(self._yoko.output_en, booleen)
    def set_range(self,val):
        if self._debug :
            print("set Yoko range : {}".format(val))
        else :
            set(self._yoko.range,val)
    """
        Extra behavior
    """
    def set_init_state(self,max_val):
        self.set(0)
        self.set_output(True) 
        self.set_range(max_val)
    def set_close_state(self):
        self.set(0)
        self.set_output(False) 
    def set_and_wait(self,V,Waittime = 0.4):
        self.set(V)
        _time.sleep(Waittime)  # Waiting until voltage is stable
        
class Guzik_wrapper(Pyhegel_wrapper):
    """
        ----------
        Last update
        
        Added dummy config
        
        ----------
        Description
        
        This is a wrapper arround pyhegel's guzik class
        It doest not fallow pyhegel's convention because I could not solve 
        some issues that arose when trying to write a child class for 
        instruments.guzik_adp7104
        
        What it does :
            - Makes sure that only one instance of instruments.guzik_adp7104 exists
            - Add some custom behavior
            
        Options :
            - The debug option bypasse the actual initialization of the card
                this is used to avoid the long ~10 sec wait time when the goal is juste to find typos and runtime bugs...
        
        You can get direct access to the pyhegel object using
            get_pyhegel_instance()
            
        Example :
            gz = Guzik_wrapper() 
            gz.config(stuff..) 
            
            new_gz = Guzik_wrapper()  # Points to the same instance of instruments.guzik_adp7104
            
            Guzik_wrapper.couter # returns [2]
            
            del new_gz  
            
            Guzik_wrapper.couter # returns [1]
            
            pyhegel_gz = gz.get_pyhegel_instance() # This can be used with pyhegel functions
            
            data    = gz.get() # Gets the data
            snippet = gz.get_snippet() # Gets the first 1000 points. Does not mess with the config
            
            del gz  # The pyhelgel object gets deleted
            
        Undefined behavior :
            - If you try to make an instance of instruments.guzik_adp7104() 
            not using this class while this class as already an existing instance/object.
            
        Knowned bugs :
            - 
        Todos :
            - Add get_or_getcache
            - Can I make _is_debug work in a way that is actually usefull (i.e. reduce repetition in the code and the clutter)?
                - Also in config() I make an instance of dummy_data but I'm not covering the case where more than one channel is declared
                - This last problem propagates to the get function since idk what is the entended return format for the data when more than one channel is declared
        Bugs :
            - 
    """
    __version__     =  {'Guzik_wrapper':0.4}
    __version__.update(Pyhegel_wrapper.__version__)
    __dummy_config__ = {'conv_resolution':_np.r_[1.0]}
    _gz_instance =   [] 
    counter     =   [0]     
    def __init__(self,debug=False,verification_data=None ):
        """
            Comment to be eventually deleted and info integrated to documentation
            
            Made the constructor more explicit for keywords so that the signature would carry more info
            
            Added verification_data 
                if None 
                    gz bahaves normally
                else : 
                    it is an array of shape (N,l_data) given by the user
                    and get will iterate through that array. So effectivelly it 
                    bypasses the get function and return a predertermined vector.
                    Additionnaly the get function is going to get a force argument to
                    force the real aquisition if verification_data is used.
        """
        self.counter[0] += 1
        self._debug  = debug
        self._v_data = verification_data
        self.__load_guzik__()
        self._gz = self._gz_instance[0] 
        self._config_inputs = dict() # used to save config calls
    def __del__(self):
        self.counter[0] -= 1 
        self._gz = None 
        if self.counter[0] <=0 :
            del self._gz_instance[0] 
            self._gz_instance.pop 
    def __load_guzik__(self):
        """
            This ensures that only one instance of instruments.guzik_adp7104 exist
            for all objects of the Guzik_wrapper class
        """
        if not self._debug : 
            try:
                if not isinstance(self._gz_instance[0], instruments.guzik.guzik_adp7104):
                    print("\nLoading guzik :")
                    self._gz_instance.append( instruments.guzik_adp7104() )
            except:
                print("\nLoading guzik :")
                self._gz_instance.append( instruments.guzik_adp7104() )
        else :
            self._gz_instance.append(None)
    def config(self, channels=None, n_S_ch=1024, bits_16=True, gain_dB=0., offset=0., equalizer_en=True, force_slower_sampling=False, ext_ref='default', _hdr_func=None):
        if channels is not None :
            self._config_inputs = dict(channels=channels, n_S_ch=n_S_ch, bits_16=bits_16, gain_dB=gain_dB, offset=offset, equalizer_en=equalizer_en, force_slower_sampling=force_slower_sampling, ext_ref=ext_ref, _hdr_func=_hdr_func)
        if not self._debug : 
            return self._gz.config(channels,n_S_ch,bits_16,gain_dB,offset,equalizer_en,force_slower_sampling,ext_ref,_hdr_func)
        else :
            self._dummy_data = _np.zeros((1,n_S_ch),dtype='int16') # this doesn't cover all the corner cases
            return Guzik_wrapper.__dummy_config__
    def read_config(self):
        if not self._debug : 
            return self._gz._read_config()
        else :
            return None
    def get_config_inputs(self):
        """
            Return a dict of the last used inputs of self.config(...)
        """
        return self._config_inputs
        
    def get_mv_per_bin(self):
        if not self._debug : 
            return _np.r_[self._gz._read_config()['conv_resolution']]*1000
        else :
            return 1.0
    def get_center_bin(self):
        """
            Subtract this amount to get the equivalent int number
        """
        if not self._debug : 
            return _np.r_[self._gz._read_config()['conv_offset']]
        else :
            return 0.
    
    def get(self,force_aquisition=False,force_verification=False):
        if ( ( (not self._debug) and (not force_verification) ) or force_aquisition ): 
            return get(self._gz)
        elif force_verification :
            print("get verification data ")
            return self._v_data
        else :
            print("get Guzik ")
            return self._dummy_data[0,:]
    def get_pyhegel_instance(self):
        return self._gz
    def get_snippet(self,snipsize=1000):
        if not self._debug : 
            return self.get()[:snipsize]
        else :
            return self._dummy_data[0,:snipsize]
    def quick_histogram_int16(self,channel=None,n_threads=32):
        hist = _Hist_int16_t(n_threads,bit=16) # Bug : not using bit=16 does work but crashes in accumulate
        if not self._debug : 
            if channel ==None:
                data = self.get()
            elif (channel>0) and (channel<5) :
                data = self.get()[channel-1]
            else:
                raise Exception(" Channel must be in 1 to 4.")
        else :
            data = self._dummy_data
        hist.accumulate(data)
        fig, axs = _plt.subplots(1,1)
        axs.plot(arange(-2.0**15,2.0**15),hist.get())
        axs.set_xlim(0,1023)
        del hist
      
            
        # Instanciate X
        X = TimeQuad_uint64_t(R,dt,l_data,kernel_conf,Filters,g,alpha,l_fft,n_threads)
        
        # Instanciate Hs
        shape = X.ks().shape[:-1]
        Hs    = build_array_of_objects( shape , _Hist_double , *(nb_of_bin,n_threads,max) ) 
        
        # Get data
        if not self._debug : 
            data = self.get()
        else :
            data = self._dummy_data
        
        X.execute( data ) # force the initialization of memory
        quads = X.quads()
        
        # Accumulate in histograms 
        quads_shape = sub_flatten ( quads, axis = -1 ) # this modifies input arr
        for quad, H in zip(quads,Hs.flat):
            H.accumulate(quad)
        # Restore shapes
        quads.shape = quads_shape
        
        # Plot all hist
        fig, axs = _plt.subplots(1,1)
        for i,H_quad in enumerate(Hs):
            for H_filters, label in zip( H_quad, Labels) :
                if i == 0 :
                    axs.plot( _np.linspace(- max,max, nb_of_bin ) ,H_filters.get() ,label = label)
                else :
                    axs.plot( _np.linspace(- max,max, nb_of_bin) ,H_filters.get() ,ls=':')
        axs.legend()
        # Restor gz
        self.config(**kw_save)
        
        # Free memory
        del quads
        del X
        del Hs

class PSG_wrapper(Pyhegel_wrapper):
    """
        Todos :
            - Implement basic behaviour
    """
    __version__ = {'PSG_wrapper':0.2}
    __version__.update(Pyhegel_wrapper.__version__)
    name_addr = {  'PSG_E8257D':'GPIB0::19::INSTR',}
    default_freq    = 20.0e9
    default_ampl    = -135.0
    default_rf_en   = False
    
    def __init__(self,visa_addr=None,nickname='PSG_E8257D',**options):
        self._set_options(**options)
        if self._debug :
            self._psg       = None
        else :
            visa_addr = visa_addr if visa_addr else PSG_wrapper.name_addr[nickname]
            self._psg = instruments.agilent_rf_PSG(visa_addr=visa_addr)
            self.set_init_state()
    def __del__(self):
        self._set_options(**options)
        if self._debug :
            pass
        else :
            self.set_close_state()
    """
        Wrapper of existing behavior
    """
    def _set_options(self,**options):
        self._debug  = options.get('debug')
    def get(self):
        if self._debug :
            return (PSG_wrapper.default_freq,PSG_wrapper.default_ampl)
        else :
            f   = get(self._psg.freq_cw)
            dbm = get(self._psg.ampl)
            return (f,dbm)
    def set_freq(self,f):
        if self._debug :
            print("set PSG : {:0.2f}[GHz]".format(f*1.e-9))
        else :
            set(self._psg.freq_cw,f)
    def set_ampl(self,dBm):
        if self._debug :
            print("set PSG : {:0.2f}[dBm]".format(dBm))
        else :
            set(self._psg.ampl,dBm)           
    def set_output(self,booleen):
        if self._debug :
            pass
        else :
            set(self._psg.rf_en, booleen)     
    """
        Extra behavior
    """
    def set_init_state(self):
        # PSG is not coherent when changing freq hence I dont change it here to stay coherent when constructor and destructor are called.
        # 
        #self.set_freq(PSG_wrapper.default_freq) 
        self.set_ampl(PSG_wrapper.default_ampl)
        self.set_output(PSG_wrapper.default_rf_en)
    def set_close_state(self):
        # PSG is not coherent when changing freq hence I dont change it here to stay coherent when constructor and destructor are called.
        #self.set_freq(PSG_wrapper.default_freq)
        self.set_ampl(PSG_wrapper.default_ampl)
        self.set_output(PSG_wrapper.default_rf_en)
        
        
