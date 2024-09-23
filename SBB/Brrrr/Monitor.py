#!/bin/python
# -*- coding: utf-8 -*-
#https://www.geeksforgeeks.org/create-a-watchdog-in-python-to-look-for-filesystem-changes/

import os as _os
import time as _time
import logging as _log                 
from watchdog.observers import Observer as _Observer
from watchdog.events import FileSystemEventHandler as _FSEH

def RunOnEvent(commands=('echo "Hello World"',),Dir=".",pull_time=1,silent=True):
    _log.basicConfig(level=_log.INFO,format='%(asctime)s - %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    event_handler = LogAndRun(commands,Dir,pull_time,silent) 
    try:
        while True:
            _time.sleep(1.0)
    except KeyboardInterrupt:
        event_handler.observer.stop()
    event_handler.observer.join()
    
import os as _os
import time as _time
import logging as _log
from watchdog.observers import Observer as _Observer
from watchdog.events import FileSystemEventHandler as _FSEH

def RunOnEvent(commands=('echo "Hello World"',), Dir=".", pull_time=1, silent=True,only_for_ext='.npz'):
    _log.basicConfig(level=_log.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    event_handler = LogAndRun(commands, Dir, pull_time, silent,only_for_ext)
    try:
        while True:
            _time.sleep(1.0)
    except KeyboardInterrupt:
        event_handler.observer.stop()
    event_handler.observer.join()

class LogAndRun(_FSEH):
    """
    Runs a list of commands when files with .npz extension are created or modified
    
    Will only trigger for file ending with only_for_ext
    """
    command = ""
    def __init__(self, commands, dir, pull_time, silent,only_for_ext):
        super(LogAndRun, self).__init__()
        self.commands = commands
        self.dir = dir
        self.monitored_dir = _os.path.abspath(self.dir)
        self.pull_time = pull_time
        self.silent = silent
        self.ext    = only_for_ext
        self.observer = _Observer()
        self.observer.schedule(self, self.dir, recursive=True)
        self.observer.start()
        print(("Monitoring directory : {}".format(self.monitored_dir)))
        if silent:
            _log.info("SILENT MODE ON. NO OUTPUT DISPLAYED.")
        
        _log.info("Waiting for event : ")

    def on_any_event(self, event):
        self.observer.unschedule_all()
        _time.sleep(self.pull_time) # time buffer so it does trigger multiple times when events are close in time
        if event.is_directory: # Do Nothing for directory changes
            pass
        elif event.event_type in ['modified', 'created'] and event.src_path.endswith(self.ext):
            print(("File {} created or modified.".format(event.src_path)))
            for cmd in self.commands:
                if self.silent:
                    cmd += " > " + _os.devnull + " 2>&1"
                _log.info("Running: " + cmd)
                _os.system(cmd)
        self.observer.schedule(self, self.dir, recursive=True)
        _log.info("Waiting for event : ")

 