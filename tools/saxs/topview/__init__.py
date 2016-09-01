import numpy as np

from .TVLog import TVLog
log = TVLog()

def namestr(name):
  assert isinstance(name, tuple), name
  s = str(name[0])
  f = name[1]
  if f is not None:
    s += "[" + str(f) + "]"
  for f in name[2:]:
    s += "::" + str(f)
  return s

class ViewError(TypeError):
  pass

class TVBase(object):
  pass

class TVData(TVBase):
  pass

_mode = "direct"

def get_mode():
  return _mode

def set_mode(mode):
  global _mode
  assert mode in ("direct", "unroll", "search_cache", "evaluate"), mode
  #print >> sys.stderr, "SET-MODE", _mode, "=>", mode
  _mode = mode

from .TVAction import TVAction, GPUAction
from .TVArray import TVArray, TVArraySeries, TVArrayList
from .TVContext import TVContext

try:
  import pycuda
  import pycuda.autoinit
except ImportError:
  pass  

from .TVCache import TVHash
