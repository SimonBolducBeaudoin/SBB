import os
s = os.path.abspath("C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin")

if os.name == "nt" and s not in os.environ["PATH"]:
  os.environ["PATH"] = s+";"+os.environ["PATH"]
  
from .acorrs_otf import *
from AutoCorr_helper import *

import Deprecated

del AutoCorr_helper