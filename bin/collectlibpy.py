import os, sys

currdir = os.path.split(__file__)[0] 
if not len(currdir): currdir = "."
so = currdir+ os.sep+"collectlib.so"
from ctypes import *
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

def collect_coor():
  x = (POINTER(c_double)).in_dll(lib, "x")
  ret = []
  start = 0
  for end in ieins[:nlig]:
    coor = []
    for n in range(start, end):
      c = x[3*n],x[3*n+1],x[3*n+2]
      coor.append(c)
    ret.append(coor)  
    start = end
  return ret

def collect_all_coor():
  x = (POINTER(c_double)).in_dll(lib, "x")
  coor = []
  for n in range(ieins[nlig-1]):
    c = x[3*n],x[3*n+1],x[3*n+2]
    coor.append(c)
  return coor

  
