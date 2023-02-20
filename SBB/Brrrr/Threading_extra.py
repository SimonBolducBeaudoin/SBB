#!/bin/python
# -*- coding: utf-8 -*-

# https://stackoverflow.com/questions/6893968/how-to-get-the-return-value-from-a-thread

import sys as _sys
from threading import Thread as _Thread

__all__ = ["ThreadWithReturnValue","parralel_bunch"]

if _sys.version_info[0] ==2 : # python 2
    class ThreadWithReturnValue(_Thread):
        """
        Join will now return a value
        """
        def __init__(self, group=None, target=None, name=None,
                     args=(), kwargs={}, Verbose=None):
            _Thread.__init__(self, group, target, name, args, kwargs, Verbose)
            self._return = None
        def run(self):
            if self._Thread__target is not None:
                self._return = self._Thread__target(*self._Thread__args,**self._Thread__kwargs)
        def join(self):
            _Thread.join(self)
            return self._return
elif _sys.version_info[0] ==3 : #python 3
    class ThreadWithReturnValue(_Thread):
        """
        Join will now return a value
        """
        def __init__(self, group=None, target=None, name=None,
                     args=(), kwargs={}, Verbose=None):
            _Thread.__init__(self, group, target, name, args, kwargs)
            self._return = None
        def run(self):
            if self._target is not None:
                self._return = self._target(*self._args,**self._kwargs)
        def join(self, *args):
            _Thread.join(self, *args)
            return self._return
else :
    raise Exception("Python version not compatible")

def parralel_bunch(funcs,SERIAL=False,args=None,kwargs = None):
    """
    Automatic parralel threads call of a list of funcs with the same args and kwargs.
    """
    if SERIAL : # for debug purposes
        if args and kwargs :
            return [ f(*args,**kwargs) for f in funcs ]
        elif args :
            return [ f(*args) for f in funcs ]
        elif kwargs :
            return [ f(**kwargs) for f in funcs ]
        else :
            return [ f() for f in funcs ]
    else:
        if args and kwargs :
            threads = [ ThreadWithReturnValue(target = f, args = a ,kwargs = kw) for f, a,kw in zip(funcs,args,kwargs) ]
        elif args :
            threads = [ ThreadWithReturnValue(target = f, args = a )             for f, a in zip(funcs,args) ]
        elif kwargs :
            threads = [ ThreadWithReturnValue(target = f, kwargs = kw)           for f, kw in zip(funcs,kwargs) ]
        else :
            threads = [ ThreadWithReturnValue(target = f ) for f in funcs ]
        for th in threads :
            th.start()
        return [ th.join() for th in threads ]