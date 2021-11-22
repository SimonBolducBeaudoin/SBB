#!/bin/env/python
#! -*- coding: utf-8 -*-

import time
import numpy
import itertools
import os

import pkg_resources  # part of setuptools
_SBB_version = pkg_resources.require("SBB")[0].version
del pkg_resources

import git
repo = git.Repo(search_parent_directories=False)
sha = repo.head.object.hexsha
del git
print "SBB version number : {}".format(_SBB_version)
print "Repository : {}".format(repo)
print "Hash : {}".format(sha)


# DEPRECATED LINE from General_tools  import compute_differential

####################
#  Utilities #
####################

def dict_to_attr(self,dict,user_key_to_attribute_key={}):
    """
        Used in a constructor to automatically set many attributes defined in a dictionnary
    """
    conv = user_key_to_attribute_key 
    for key in dict :
        setattr(self,conv[key],dict[key]) if key in conv else setattr(self,key,dict[key])

####################
# Pyhegel Utilities #
####################
class timer(object):
    def __init__(self,size=1):
        self.timer = numpy.zeros((2,size))
    def watch(self):
        return time.time()  
    def tic(self,index=0):
        self.timer[0,index] = self.watch()
    def toc(self,index=0):
        self.timer[1,index] = self.watch()-self.timer[0,index]
    def durations(self):
        return self.timer[1,:]
        
class logger(object):
    """
        Times an experiment 
        and print events and progress
        begin
            loop
                loop_events
        end
        
        log_dict (example)
        {
            'open'          : 'Lauching experiment'
            'close'         : 'Experience over : \n \t Estimated time [s] {:04.2F} \n \t Real time {:04.2F} [s]'
            'loop_sizes'    : (1st_loop,2nd_loop,...)
            'loop'          : 'Progess {:03}/{:03d} \t last_loop time: {:04.2F} [s] '
            'events'        : (None by default)
                            : (format example) # The first parameter is an index for the event since dictionnary are now guarentied to preserve order
                                {
                                    "Event_name" : "format_ready_str"
                                }
                            : (example) 
                                {
                                "Aquisition": "Acquisition : {:04.2F} [s]", 
                                "Computing" : "Computing : {:04.2F} [s] "
                                }
            'rate'          : ( l_data*1.0e-9.0 ,"Rate : {:04.2F} [GSa/s] " ),
        }
        times_est   = (experiment_time,loop_time)
        
        Todos :
            - Modifiy the _print() to save (or not) the string of log into a file
          Bugs :
            -
    """   
 
    _default_log = \
    {
        'loop_sizes': (1,),
        'open'      : 'Lauching experiment' ,
        'close'     : 'Experience over : \n \t Estimated time [s] {:04.2F} \n \t Real time {:04.2F} [s]' ,
        'loop'      : 'Progess {:03}/{:03d} \t last_loop time: {:04.2F} [s]' ,
        'conditions': '{:04.1F}'
    }
    _allowed__user_key = ('loop_sizes','open','close','loop','rate','events')
    _user_key_to_attribute_key = \
    {
        'loop_sizes'    : '_loop_sizes'  ,
        'open'          : '_open_str'    ,
        'close'         : '_close_frmt'  ,
        'loop'          : '_loop_frmt'   ,
        'is_rate'       : '_is_rate'     ,
        'rate'          : '_event_rate'  ,
        'is_event'      : '_is_event'    ,
        'events'        : '_events'      ,
        'conditions'    : '_conditions_frmt'
    }
    _default_time_est = (1.0,)
    def __init__( self, time_estimates , log_dict , **options):
        self._attributes_instanciation(time_estimates,log_dict)
        self.logger_indent      = options['indent'] if options.has_key('indent') else 0
    def indent(self,n):
        self.logger_indent = n
    def _print(self,s):
        ss =''
        for i in range(self.logger_indent) :
            ss += '\t'
        print ss + s
    def _attributes_instanciation(self,time_estimates , log_dict):
        if log_dict == () or log_dict == None : # Set the default log
            dict_to_attr(self,self._default_log,self._user_key_to_attribute_key)
        else :
            self._log_dict_to_attributes(log_dict,self._user_key_to_attribute_key)
        self._event_attributes_instanciation(log_dict)
        self._set_time_estimate(time_estimates)
    def _log_dict_to_attributes(self,log_dict,conv):
        for key in self._default_log :  # set the attributes of the given log_dict or the default one when it doens't exist
            setattr(self,conv[key],log_dict[key]) if key in log_dict else setattr(self,conv[key],self._default_log[key])
    def _set_time_estimate(self,time_estimates):
        self.time_estimate =  self._default_time_est if  time_estimates  == () or time_estimates == None else time_estimates
    def _event_attributes_instanciation(self,log_dict):
        self._is_rate    =   'rate' in log_dict
        if self._is_rate :
            self._event_rate = log_dict['rate']
        self._experiment = timer()
        self._experiment.tic()
        self._experiment.toc()
        self._loop       = timer(len(self._loop_sizes))
        self._loop.tic()
        self._loop.toc()
        self._is_event   = 'events' in log_dict
        if self._is_event :
            self._events_dict = log_dict['events']
            self._events_len  = len(self._events_dict)
            self._events      = timer(self._events_len)
            self._events.tic(0)
            self._events.toc(0)
            self._conditions_frmt = log_dict['conditions_frmt'] if log_dict.has_key('conditions_frmt') else ('{: .1f}',)  
    def open(self):
        self._experiment.tic()
        for i in range(len(self._loop_sizes)-1):
            self._loop.tic(i)
        if self._is_event :
            self._events.tic(0)
        self._print(self._open_str)
    def close(self):
        self._experiment.toc()
        if self._is_event :
             self._events.toc(self._events_len-1)
        total_t = self._experiment.durations()[0]
        self._print((self._close_frmt).format(self.time_estimate[0],total_t))
    def loop(self,loop_index,loop_icrmnt):
        self._loop.toc(loop_index)
        loop_time  = self._loop.durations()[loop_index]
        self._loop.tic(loop_index)
        self._print (self._loop_frmt.format(loop_icrmnt, self._loop_sizes[loop_index] , loop_time))
    def event(self,index):
        self._events.toc(index)
        self._events.tic((index+1)%self._events_len)  
    def _build_events_frmt(self):
        drtns = self._events.durations()
        s = ''
        tab = ' | '
        d = self._events_dict
        for i,e in enumerate(d):
            s += tab
            s += e.format(drtns[i])
        if self._is_rate:
            s += tab
            num = self._event_rate[0]
            frmt = self._event_rate[1]
            dur = drtns[-1]+drtns[0]
            if dur <= 0 :
                s += 'Nan'
            else :
                s += frmt.format(num/dur)
        return s     
    def events_print(self,condition_tuple):
        s = self._build_events_frmt()
        s_tuple = '('
        for i,t in enumerate(condition_tuple) :
            try :
                cnd_frmt = self._conditions_frmt[i]
            except IndexError : 
                cnd_frmt = self._conditions_frmt[0]
            s_tuple += cnd_frmt.format(t) + ','
        s_tuple += ')'
        self._print(s_tuple + '\t' + s)
    
class ExperimentErrors(Exception):
    """
        Base class for pyHegel tools errors
    """
    pass

class ExperimentWarnings(Warning):
    """
        Base class for pyHegel tools warnings
    """
    @classmethod
    def warn(cls,*arg,**kwargs):
        self = cls(*arg,**kwargs)
        s = str(self)
        name = cls.__name__
        print name + ' : \n' + s
        return self

class LoadingWarning(ExperimentWarnings):
    """
       When loading from saved data
    """
    pass
    
class VersionsError(ExperimentErrors):
    """
        Deprecated
    """
    def __init__(self,versions_current,versions_loaded):
        self.versions_current  = versions_current
        self.version_loaded = versions_loaded
        s_load     = "Loaded versions  : \n\t {}".format(versions_loaded)
        s_current  = "Current versions : \n\t  {}".format(versions_current)
        self.message    = s_load+'\n'+s_current
        super(VersionsError,self).__init__(self.message)
  
class VersionsWarning(ExperimentWarnings):
    def __init__(self,versions_current,versions_loaded):
        self.versions_current  = versions_current
        self.version_loaded = versions_loaded
        s_load     = "Loaded versions  : \n\t {}".format(versions_loaded)
        s_current  = "Current versions : \n\t  {}".format(versions_current)
        self.message    = s_load+'\n'+s_current
        super(VersionsWarning,self).__init__(self.message)

class ConditionsError(ExperimentErrors):
    """
        When the given conditions are not of the right type or shape 
        Or error regarding Conditions_logic class
    """
    
class UndefinedStateError(ExperimentErrors):
    """
       When the experiment object ends in an undefined state 
    """
    pass 
   
#################
# Pyhegel tools #
#################

class Info(object):
    """
        Contains all the info about a given experiement
            - Sweeping conditions are contained in a tuple of 1D array. 
                We assume that the experiement will sweep over the outerproduct of those array, though that behavior might be modified by the user in the experiment class.
            - n_measures represents the number of time the conditions are repeated in the experiement
            - Meta_info is a dictionnary contaning meta information about the experiement (anything really)
                A repetition keys is added to meta_info in case the experiment is repeated by the user (see: Experiment class)
        
        Also set a few defaults attributes likes self._verboe and self._test ...
    """
    __version__     = { 'Info'  : 1.1 }
    def _set_options(self,options):
        self._options                   = options
        self._verbose                   = options.get('verbose')
        self._test                      = options['test']    if options.has_key('test') else True 
    def _set_conditions(self,conditions):
        self._conditions                = conditions
        if type(conditions[0]) != int :
            raise ConditionsError('n_measures should be int')
        self._n_measures                = conditions[0]         # The first element of the tuple is the number of repetions
        self._conditions_core_loop_raw  = conditions[1:]        # The 2nd   element ...          is an list of 1D array    
        self.conditions                 = self.get_conditions() # Public copy of the experimental conditions
        self.n_measures                 = self.get_n_measures() # Public copy of _n_measure
    def _set_meta_info(self,meta_info):
        self._meta_info                 = meta_info
        self._meta_info['repetitions']  = meta_info['repetitions'] if meta_info.has_key('repetitions') else 0
    def _build_attributes(self):
        pass
    def get_conditions(self):
        """
        Returns a tuple of 1D array corresponding to the experimental conditions
        """
        return self._conditions_core_loop_raw
    def get_n_measures(self):
        return self._n_measures

class Accretion(object):
    """
        TLDR :
        Todos :
            - __doc__
        Bugs :
        Knowed issues :
            Using numpy savez to save data makes it such that allow_pickle option must be used.
            Dictionnary therefore must be loaded in a particular way and list and tuples are converted
            to np_array.
    """
    __version__     = { 'Accretion'  : 1.0 }
    @classmethod
    def description(cls):
        print cls.__doc__
    def __init__(self, exp):
        self.exp = exp
    #############
    # User interface #
    #############
    # def update_analysis   : (below)
    # def save_data         : (below)
    # def load_data         : (below)
    # def get_data_dict     : (below)
    #############
    # Utilities #
    #############
    def get_data_dict(self):
        return self._data
    ############
    # Reduction/Analysis #
    ############
    def _compute_analysis(self):
        pass
    def _build_data(self):
        return {} # dummy default behavior
    def _update_data(self):
        self._data = self._build_data()
    def update_analysis(self,**kwargs):
        self._compute_analysis(**kwargs)
        self._update_data()
    #############
    # Save/load #
    #############
    def save_data(self,path_save,prefix='acc_'):    
        time_stamp                  = time.strftime('%y%m%d-%H%M%S') # File name will correspond to when the experiment ended
        filename                    = prefix+'{}.npz'.format(time_stamp)
        to_save                     = self._data
        to_save['_versions_saved']  = self.__version__
        numpy.savez_compressed(os.path.join(path_save,filename),**to_save)
        print "Data saved \n \t folder : {} \n \t {}".format(path_save,filename) 
    def _load_data_dict(self,data_dict):
        dict_to_attr(self,data_dict)
        self._data  = data_dict
    @classmethod
    def load_data(cls,folder,filename):
        """
            To load data create an experiment object using
            experiment = class_name().load_data(folder,filename)
        """
        data                    = numpy.load(folder+'\\'+filename,allow_pickle=True)
        data                    = dict(data)
        self._load_data_dict(data_dict)
        return self
        
class Analysis(Info):
    """
        TLDR :
            - A class for the analysis of experiments (see: class Experiment)
            - can be constructed either 
                - from an experiment class object
                - or from conditions,meta_info,data,**options
                - or from saved data
            - It will have all the attibutes defined by conditions,meta_info,data and options
            - It will check for version compatibility
            - User will write a class that inherits from this one for specific experiment
        Todos :
            - Update analyse in init ...
            - Change update analysis to update
        Bugs :
        Knowed issues :
            Using numpy savez to save data makes it such that allow_pickle option must be used.
            Dictionnary therefore must be loaded in a particular way and list and tuples are converted
            to np_array.
    """
    __version__     = { 'Analysis'  : 1.0 }
    @classmethod
    def description(cls):
        print cls.__doc__
    def __init__(self,conditions,data_dict,meta_info=None,**options):
        """
            - conditions    : example == (n_measures, nd_array_core_loop_cnds) 
            - devices       : example == (guzik,yoko,)
            - meta_info     : example == (l_kernel,l_data,R_jct,dt,)
            - option        : kwarg depends on the specific child class
        """   
        self._set_options(options)
        self._set_conditions(conditions)
        self._set_meta_info(meta_info)
        self._build_attributes()
        self._load_data_dict(data_dict)
    #############
    # User interface #
    #############
    # def update_analysis   : (below)
    # def save_data         : (below)
    # def load_data         : (below)
    # def get_data_dict     : (below)
    ######################
    # Analysis Utilities #
    ######################
    @staticmethod
    def SE(mu2k,muk,n):
        """ 
            Voir notes Virally Central limit theorem
            Computation of the standard error for the moment of order K
            mu2k : is the moment of order 2 k
            muk  : is the moment of order k
            If these moments are not centered then the definition is good for none centered moment
            Idem for centered moment
        """
        return numpy.sqrt(numpy.abs(mu2k-muk**2)/float(n))
    #############
    # Utilities #
    #############
    def _compute_n_measure(self):
        return 1 if self._test else self._n_measures
    def _n_measure_total(self):
        return self._n_measures*(1+self._meta_info['repetitions'])
    def get_data_dict(self):
        return self._data
    ############
    # Reduction/Analysis #
    ############
    def _compute_analysis(self):
        pass
    def _build_data(self):
        return {} # dummy default behavior
    def _update_data(self):
        self._data = self._build_data()
    def update_analysis(self,**kwargs):
        self._compute_analysis(**kwargs)
        self._update_data()
    #############
    # Save/load #
    #############
    def save_data(self,path_save,prefix='anal_'):    
        time_stamp                  = time.strftime('%y%m%d-%H%M%S') # File name will correspond to when the experiment ended
        filename                    = prefix+'{}.npz'.format(time_stamp)
        to_save                     = self._data
        to_save['_versions_saved']  = self.__version__
        to_save['_options']         = self._options
        to_save['_conditions']      = self._conditions
        to_save['_meta_info']       = self._meta_info
        numpy.savez_compressed(os.path.join(path_save,filename),**to_save)
        print "Data saved \n \t folder : {} \n \t {}".format(path_save,filename) 
    def _load_data_dict(self,data_dict):
        dict_to_attr(self,data_dict)
        self._data  = data_dict
    def _check_cls_vs_data_versions(self):
        try : 
            versions_saved = self._versions_saved
        except AttributeError :
            versions_saved = None
        version = self.__version__
        try :
            version.pop(type(self).__name__)
            version.pop(Analysis.__name__)
        except KeyError :
            LoadingWarning.warn('Analysis class did not contain its own version number.')
        if ( version != versions_saved ) and versions_saved:
            VersionsWarning.warn(version,versions_saved)
    @classmethod
    def load_data(cls,folder,filename):
        """
            To load data create an experiment object using
            experiment = class_name().load_data(folder,filename)
        """
        data                    = numpy.load(folder+'\\'+filename,allow_pickle=True)
        data                    = dict(data)
        try :
            conditions          = data.pop('_conditions')
        except KeyError as er:
            raise LoadingWarning('Missing key '+er.args[0]+' in loaded data.')
        try :
            meta_info           = data.pop('_meta_info')[()]        # to load a dict saved by numpy.savez
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            meta_info = None
        try :
            options             = data.pop('_options')[()]          # to load a dict saved by numpy.savez
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            options = None
        self                    = cls(conditions,data,meta_info,**options)
        self._check_cls_vs_data_versions()
        try :
            self._versions_saved= data.pop('_versions_saved')[()] 
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            self._versions_saved = None
        return self
    @classmethod 
    def load_exp(cls,exp):
        conditions  = exp._conditions
        meta_info   = exp._meta_info
        options     = exp._options
        data        = exp.get_data_dict()
        return cls(conditions,data,meta_info,**options)
    
class Experiment(Analysis):
    """
        TLDR :
            (Pseudo code (not up to date): Experiment) 
                experimentals condition, meta_info, pyhegel devices and options are given to the constructor
                
                child = children_class(conditions,devices,meta_info,op1=True,op2=True,op3=x,)
                
                The constructor initialize internal variable to a safe states.
                It instanciate computing object (autocorrelation,convolutions,...)
            
                child.measure() calls the measurment loop :
                    all_loop_open
                    for n  in main_iterator :               # controls repetitions using the default main_iterator (if not overwritten)
                        repetition_loop_start(n)            # code that as to be executed before each measurement_loop
                        for ... core_loop_iterator :        # using the user defined core_loop_iterator
                            loop_core(i_tuple,cndtns_tuple) # note that i and cndtns are tuples
                        repetition_loop_end(n)
                    all_loop_close()
                
                child.update_analysis()                     # compute the relevant physical data
                
                child.fig...()                              # standard plot(s) for this experiment
                
                child.save_data(folder)                     # saves data_timestamps.npz into folder
                                                                version number are saved with the data
                                                                
                child.load_data(folder,filename)            # construct a child object from saved data
                
                data = child.get_data()                     # return child._data a list of array defined by the user
                                                                # The user can modify that data and update_analysis() followed by plots to tinker with the data
        
        - New with version 1.0
            - There is now a Analysis dedicated class
            
        HOW TO WRITE THE CHILD CLASS (Not up to date):
        This class is written to give flexibility and reduce the amount of code that the user as to writte for the child class.
        This section describes what the user as to write in the child class
        - __init__ function
            You should not have to write a __init__ function for the child class
            The behavior of __init__ is modified by overwriting the following methods
            - set functions
                - _set_options       : setting attributes that modifies behaviors of other methods
                - _set_conditions    : sets self.conditions (mendatory) 
                                        and unrolls 
                                        self._n_measures = conditions[0] (mendatory)
                - _set_meta_info     : sets self._meta_info (mendatory) and unrolls it
                - _build_attributes  : build other attributes (Falcultative)
                - _set_devices       : sets self.devices (mendatory) 
                                        , unrolls it
                                        and sets each devices into their initial sate (Good pratice)
                - _init_objects      : initialize memory for large objects that will be used recurently (Falcultative)
                - _build_data        : helps to build self._data dict
                - _init_log          : defines the dehavior of the timer (Falcultative)
        , defining them is in a sense falcultative as the child class will still work if some or all of those are not defined 
        but it may lead to undefined bahavior
        - __del__ function
            It is not mendatory to write a __del__ function but a good pratice is 
            to set all devices to a safe/stable close state using __del__
        - Utilities
            Contains methods that are usefull all-around
        - loop behavior
            Calling the Experiment.measure() method will launch a measurment loop 
            
            (Pseudo code : measurment loop)
                for n  in main_iterator :               # controls repetitions using the default main_iterator (if not overwritten)
                    loop_before()                      # code that as to be executed before each measurement_loop
                    for ... core_loop_iterator :        # using the user defined core_loop_iterator
                        loop_core(i_tuple,cndtns_tuple) # note that i and cndtns are tuples
                loop_close()                            # code that is executed after all the repetitions are over
            
            see : Experiment._repetition_loop_iterator for the exact implementation
            
            User minimally defines in child class :
                - child._core_loop_iterator(self)
                - child._loop_core(self,index_tuple,condition_tuple)
            But can also define
                - child._repetition_loop_iterator(self)
                - child._repetition_loop_start(self,n)
                - child._repetition_loop_end(self,n)
                - child._all_loop_close(self)
                and _log events (to change the log behavior) 
                                
            The _super_enumerate utility is meant to help writing _core_loop_iterator      
                        
        - analysis 
            User can call child.update_analysis() to analyse data after all the repetitions are completed
            For this to work user as to define
                - _compute_reduction()
                    Is converting the data contained in the computing objects to esally saveable np.arrays.
                - _compute_analysis()
                    Does the rest of the analysis job (i.e. calling function to convert np.arrays to np.arrays)
                - _update_data()
                    builds self.data list (containts its structure) from internal variables/attributes
            Writting those function defines the proper behavior when 
            updating the analysis and loading from existing data and 
                
        - Plot and figs
            - The plot section containt default plot that the user is probably going to want
            - every fig function must return a fig object for outside modification
        
        - Methods and variables
            - __variables/__methods are not modifiable
            - _CAPS_.. are overwritable by the child function (used for grand-mother, mother, child situations)
        
        Todos :
            - Add a way to display parents descriptions using self.description
            - Update __doc__
            - Add reset_obj ?
                - done for some sub classes but idk if im done implementing this...
        Bugs :
            - After constructing from data structure the destructor stills try to clean devices ?
    """
    __version__     = { 'Experiment'  : 1.0 }
    __version__.update(logger.__version__)
    def __init__(self,conditions,devices,meta_info=None,**options):
        """
            - conditions    : example == (n_measures, nd_array_core_loop_cnds) 
            - devices       : example == (guzik,yoko,)
            - meta_info     : example == (l_kernel,l_data,R_jct,dt,)
            - option        : kwarg depends on the specific child class
        """   
        self._set_options(options)
        self._set_conditions(conditions)
        self._set_meta_info(meta_info)
        self._build_attributes()
        self._SET_devices(devices)
        self._INIT_objects()
        self._INIT_log()
    def _set_options(self,options):
        super(Experiment,self)._set_options(options)
        self._data_from_experiment  = not(options['loading_data']) if options.has_key('loading_data') else True
        self._estimated_time_per_loop = options['estimated_time_per_loop'] if options.has_key('estimated_time_per_loop') else 1.0
        self._debug                 = options.get('debug')     
    def _SET_devices(self,devices):
        if devices == None or devices == () or self._data_from_experiment == False :
            self._devices = None
        else :
            self._devices  =   devices
            self._set_devices(devices)
    def _set_devices(self,devices):
        pass
    def _INIT_objects(self):
        if self._data_from_experiment == False :
            pass
        else :
            self._init_objects()
    def _init_objects(self):
        pass
    def _INIT_log(self):
        if self._data_from_experiment == False :
            pass
        else :
            self._init_log()  
    def _init_log(self):
        self._log            =   logger((),()) # Default timer    
    #############
    # User interface #
    #############
    def measure(self):
        self._measurement_loop()
    def repeat_measure(self,n_repetitions = 1):
        """
            In the current state of the class if I do
            n_measures = n_measures + condition[0]
            the repetition is going to be longer than the previous execution
            to keep track of the number of total number of repetition im going to temporarly use 
            a new variable, but this means that conputations using n_measures (like all _std ) wont be correct.
            
            Repetitions are saved/loaded automatically in meta_info
            Todos :
                - Update the std calculations to include repetitions
                - Anything else regarding that feature ?
            Bugs :
        """
        for  rep in range(n_repetitions):
            self._meta_info['repetitions'] += 1
            self._measurement_loop()
    # def update_analysis   : (below)
    # def save_data         : (below)
    # def load_data         : (below)
    #############
    # Utilities #
    #############
    @staticmethod
    def _super_enumerate(*args):
        """
            Args are all 1D array
        """
        if len(args) == 0 : # called empty
            return iter(()) , iter(())
        index_vec = ()
        for a in args :
            index_vec += ( range(len(a)) , )
        return itertools.product(*index_vec) , itertools.product(*args)
    def _set_all_devices(self,conditions):
        """
            Should always be implemented
        """
        pass
    def _set_devices_to_close_state(self):
        """
            Should always be implemented
        """
        if self._verbose :
            print '_set_devices_to_close_state'
        for dev in self._devices :
            dev.set_close_state()
    #################
    # Loop behavior #
    #################
    # all_loop_open
    # for n  in main_iterator :               
        # repetition_loop_start(n)            
        # for ... core_loop_iterator :        
            # loop_core(i_tuple,cndtns_tuple) 
        # repetition_loop_end(n)
    # all_loop_close()
    def _all_loop_open(self):
        self._log.open()
    def _repetition_loop_iterator(self):
        return range(self._compute_n_measure())
    def _repetition_loop_start(self,n):
        self._log.loop(0,n) 
    def _core_loop_iterator(self):
        return Experiment._super_enumerate(*self._conditions_core_loop_raw) # by default no modification on the raw input
    def _loop_core(self,index_tuple,condition_tuple):
        self._log.events_print(condition_tuple)
    def _repetition_loop_end(self,n):
        pass
    def _all_loop_close(self):
        if self._data_from_experiment :
            if self._verbose :
                print 'Setting device to close state'
            self._set_devices_to_close_state()
        self._log.close()
    def _measurement_loop(self):
        main_it = self._repetition_loop_iterator()
        self._all_loop_open()
        for n  in main_it :
            index_it , condition_it = self._core_loop_iterator()
            self._repetition_loop_start(n)
            for index_tuple, condition_tuple in zip(index_it,condition_it):
                self._loop_core(index_tuple,condition_tuple)
            self._repetition_loop_end(n)
        self._all_loop_close()
    ############
    # Reduction/Analysis #
    ############
    def _compute_reduction(self):
        pass
    def update_analysis(self,**kwargs):
        return self._update_analysis_from_aquisition(**kwargs) if self._data_from_experiment else self._update_analysis_from_load(**kwargs)  
    def _update_analysis_from_aquisition(self,**kwargs) :
        self._compute_reduction()
        self._compute_analysis(**kwargs)
        self._update_data()
    def _update_analysis_from_load(self,**kwargs):
        self._compute_analysis(**kwargs)
        self._update_data()
    #############
    # Save/load #
    ############# 
    def save_data(self,path_save,prefix='exp_'):    
        super(Experiment,self).save_data(path_save,prefix=prefix)
    def _check_cls_vs_data_versions(self):
        try : 
            versions_saved = self._versions_saved
        except AttributeError :
            versions_saved = None
        version = self.__version__
        if not ( version == versions_saved ):
            VersionsWarning.warn(version,versions_saved)
    @classmethod
    def load_data(cls,folder,filename):
        """
            To load data create an experiment object using
            experiment = class_name().load_data(folder,filename)
        """
        data                    = numpy.load(folder+'\\'+filename,allow_pickle=True)
        data                    = dict(data)
        try :
            conditions          = data.pop('_conditions')
        except KeyError as er:
            raise LoadingWarning('Missing key '+er.args[0]+' in loaded data.')
        devices                 = None
        try :
            meta_info           = data.pop('_meta_info')[()]        # to load a dict saved by numpy.savez
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            meta_info = None
        try :
            options             = data.pop('_options')[()]          # to load a dict saved by numpy.savez
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            options = None
        options['loading_data'] = True
        self                    = cls(conditions,devices,meta_info,**options)
        try :
            self._versions_saved= data.pop('_versions_saved')[()] 
        except KeyError as er :
            LoadingWarning.warn('Missing key '+er.args[0]+' in loaded data.')
            self._versions_saved = None
        self._load_data_dict(data)
        self._check_cls_vs_data_versions()
        return self          
    
class Lagging_computation(Experiment):
    """
        Initiate the conditions iterator before the core loop, but not the index iterator.
        The first conditions is set before the core loop 
        Todos : 
            - 
        Bugs :
    """
    __version__     = { 'Lagging_computation'  : 0.4 }
    __version__.update(Experiment.__version__)
    #################
    # Loop behavior #
    #################
    def _set_and_wait_all_devices(self,conditions):
        """
            Used to make sure that devices have reached stady state before
            calling the lagging core loop
        """
        pass
    def _repetition_loop_start(self,n,condition_it):
        super(Lagging_computation,self)._repetition_loop_start(n)
        self._first_conditions = condition_it.next()
        self._log.events_print(self._first_conditions)
        self._set_and_wait_all_devices(self._first_conditions)
    def _last_loop_core_iteration(self):
        pass
    def _measurement_loop(self):
        main_it = self._repetition_loop_iterator()
        self._all_loop_open()
        for n  in main_it :
            index_it , condition_it = self._core_loop_iterator()
            self._repetition_loop_start(n,condition_it)
            core_it     = iter(zip(index_it,condition_it))
            while True :
                try :
                    index_tuple, condition_tuple = core_it.next()
                    self._loop_core(index_tuple,condition_tuple)
                except StopIteration :
                    self._last_loop_core_iteration()
                    break
            self._repetition_loop_end(n)
        self._all_loop_close()

class Cross_Patern_Lagging_computation(Lagging_computation):
    """
     LAST MODIFICATION :
            self._loop_core(...)
            Modified to take facultative arguments index_it and condition_it the iterators of the core_loop
            They can be used to modify the core_loop behaviour depending on some interal values of  index_it and condition_it
    """
    class cross_enumerate(object):
        """
        An enumerator of the indexes and an iterator goes over the arguments one by one fixing all the other arguments to theirfirst value
            Example :
                cross_enumerate([1,2,3],[4,5,6]) will return consecutively
                    (1,4)
                    (2,4)
                    (3,4)
                    (1,4)
                    (1,5)
                    (1,6)
        """
        def __init__(self,*args):
            self.cndtns          = args
            try :
                self.dim             = len(args)
            except :
                raise Exception("Bad initialization")
            self.next_idx        = 0
            self.next_dim        = 0
            self.idx             = self.next_idx
            self.current_dim     = self.next_dim
        def __iter__(self):
            return self 
        def gen_idxs(self):
            index_vec = tuple()
            for a in self.cndtns :
                index_vec += ( range(len(a)) , )
            return index_vec    
        def __next__(self):
            self.current_dim    = self.next_dim
            try :
                len_current_dim     = len(self.cndtns[self.current_dim]) # wont work is empty tuple()
            except IndexError:
                raise StopIteration
            self.idx            = self.next_idx
            if self.idx == len_current_dim-1 :   # last element of current dim
                self.next_idx  = 0
                self.next_dim += 1
            else :
                self.next_idx += 1
     
            cndtn = tuple()
            for i,c in enumerate(self.cndtns):
                if i == self.current_dim :
                    cndtn +=(c[self.idx],)
                else :
                    cndtn +=(c[0],)
            return cndtn
        def next(self):
            """
                For compatibility with python 2 iterators
                see : https://stackoverflow.com/questions/5982817/problems-using-next-method-in-python
            """
            return self.__next__()
    
    # def _loop_core(self,index_tuple,condition_tuple,index_it,condition_it):
        # """
        # Loop core is defined in child class
        # For this class index_it and condition_it might be usefull
        # Copy paste this method into child class definition raplacing this comment for core_loop behaviour
        # and leaving commands bellow for log compatibility.
        # """
        # super(Croos_Patern_Lagging_computation,self)._loop_core(index_tuple,condition_tuple)
    
    class core_iterator(object):
        """
        Behaves like zip(idx_it,cdn_it) but is an iterator instead of a list
        """
        def __init__(self,idx_it , cdn_it):
            self.idx_it = idx_it
            self.cdn_it = cdn_it
        def __iter__(self):
            return self    
        def __next__(self):
            try :
                return self.idx_it.next() , self.cdn_it.next()
            except StopIteration :
                raise StopIteration
        def next(self):
            """
                For compatibility with python 2 iterators
                see : https://stackoverflow.com/questions/5982817/problems-using-next-method-in-python
            """
            return self.__next__()
    
    def _core_loop_iterator(self):
        it     = self.cross_enumerate(*self._conditions_core_loop_raw)
        it_idx = self.cross_enumerate(*it.gen_idxs())
        return it_idx, it
    
    def _measurement_loop(self):
        main_it = self._repetition_loop_iterator()
        self._all_loop_open()
        for n  in main_it :
            index_it , condition_it = self._core_loop_iterator()
            self._repetition_loop_start(n,condition_it)
            core_it = self.core_iterator(index_it , condition_it)
            for (index_tuple, condition_tuple) in core_it :
                self._loop_core(index_tuple,condition_tuple,index_it,condition_it,n)                                       
            self._last_loop_core_iteration(n)        
            self._repetition_loop_end(n)     # This feals like a repetition of self._last_loop_core_iteration is it usefull ?
        self._all_loop_close()

##################
# DEPRECATED BELOW 
##################
