import os as _os
s = _os.path.abspath("C:/cygwin64/usr/x86_64-w64-mingw32/sys-root/mingw/bin")

if _os.name == "nt" and s not in _os.environ["PATH"]:
  _os.environ["PATH"] = s+";"+_os.environ["PATH"]
  
__all__ = ["histograms","histograms_helper","moments_cumulants"]

del s
