#!/bin/env/python
#! -*- coding: utf-8 -*-

import os
s = os.path.abspath("C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin")

if os.name == "nt" and s not in os.environ["PATH"]:
  os.environ["PATH"] = s+";"+os.environ["PATH"]
  
from SBB.AutoCorr import acorrs_otf
from SBB.AutoCorr import AutoCorr_helper
from SBB.AutoCorr import Deprecated

# REMOVING UNDESIRED NAME FROM NAME SPACE
del s
del os