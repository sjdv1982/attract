import os, sys

currdir = os.path.split(__file__)[0] 
if not len(currdir): currdir = "."
so = currdir+ os.sep+"collectlib.so"
from ctypes import *
import numpy
lib = None
ieins = None

def reload():
  global lib
  if lib is not None:
    handle = lib._handle # obtain the SO handle 
    cdll.LoadLibrary("libdl.so").dlclose(handle)    
  lib = CDLL(so)
  lib.collect_init.argtypes = [c_int, POINTER(c_char_p)]

lib = CDLL(so)
lib.collect_init.argtypes = [c_int, POINTER(c_char_p)]

def collect_init(args):
  global nlig, ieins
  args =  ["collect"] + args
  args = (c_char_p * len(args))(*args)
  lib.collect_init(len(args), args )
  nlig = c_int.in_dll(lib, "nlig").value
  #maxlig = c_int.in_dll(lib, "maxlig").value
  ieins = (POINTER(c_int)).in_dll(lib, "ieins")

def check_sizes(sizes, filenames):
  assert ieins is not None
  start = 0
  for inr,i in enumerate(ieins[:len(sizes)]):
    collectsize = i-start
    if collectsize != sizes[inr]:
      raise Exception(
  "Parsing difference between collect and Python: PDB %s: %d vs %d atoms" % (filenames[inr], collectsize, sizes[inr])
  )
    start = i


def collect_iattract(args):
  args =  ["collect"] + args
  args = (c_char_p * len(args))(*args)
  lib.collect_iattract(len(args), args )
  nlig = c_int.in_dll(lib, "nlig").value
  #maxlig = c_int.in_dll(lib, "maxlig").value
  ieins = (POINTER(c_int)).in_dll(lib, "ieins")
  
def collect_next():
  ret = lib.collect_next()
  return ret

def collect_coor(copy=False):
  x = (POINTER(c_double)).in_dll(lib, "x")  
  coor = numpy.ctypeslib.as_array(x, shape=(ieins[nlig-1],3))
  if copy: coor = coor.copy()
  return coor

def collect_coor_raw(copy=False):
  x = (POINTER(c_double)).in_dll(lib, "x")  
  coor = numpy.ctypeslib.as_array(x, shape=(ieins[nlig-1]*3,))
  if copy: coor = coor.copy()
  return coor

collect_all_coor = collect_coor

  
