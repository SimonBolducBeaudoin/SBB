import os
s = os.path.abspath("C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin")

if os.name == "nt" and s not in os.environ["PATH"]:
  os.environ["PATH"] = s+";"+os.environ["PATH"]
  
from .time_quadratures import * 
from .TimeQuadrature_helper import*
import Deprecated

# REMOVING UNDESIRED NAME FROM NAME SPACE
del s
del os
del time_quadratures
del TimeQuadrature_helper
