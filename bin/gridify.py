import os
currdir = os.path.split(__file__)[0] 
if not len(currdir): currdir = "."
so = currdir+ os.sep+"_gridify.so"
from ctypes import *
lib = CDLL(so)
import numpy
lib.gridify.argtypes = [
  numpy.ctypeslib.ndpointer(dtype=numpy.float64,ndim=2,flags=("ALIGNED", "C_CONTIGUOUS")),
  c_int,
  numpy.ctypeslib.ndpointer(dtype=numpy.float64,ndim=3,flags=("ALIGNED", "F_CONTIGUOUS", "WRITEABLE")),
  c_int,
  c_int,
  c_int,
] 

def gridify(coor, grid, gridspacing, grid_origin):  
  assert isinstance(coor, numpy.ndarray)
  assert len(coor.shape) == 2 and coor.shape[-1] == 3
  grid_origin = numpy.array(grid_origin)
  assert len(grid_origin.shape) == 1 and grid_origin.shape[-1] == 3
  coor2 = numpy.array(coor, order="C")  
  coor2 -= grid_origin
  coor2 /= gridspacing
    
  dimx, dimy, dimz = grid.shape 
  lib.gridify(coor2, len(coor), grid, dimx, dimy, dimz)
  