from topview import GPUAction, TVArray
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os
from jinja2 import Template

sourcefile = os.path.split(__file__)[0] + os.sep + "/atomdensitymask.cu"

class atomdensitymask1(GPUAction):
  THREADS_PER_BLOCK = 256
  main = "atomdensitymask1"
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert accumulate
    assert len(output) == 1
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 4
    assert len(output) == 1
    d = OrderedDict(zip(("coors", "weights", "origin", "gridspacing"),input_args))
    grids = output[0]
    assert grids.shape is not None and len(grids.shape) == 4 #grids must be pre-shaped
    d["size_template"] = d["coors"].shape[1]
    d["grids"] = grids
    return d
  def get_source(self, args):    
    coors, weights, origin, gridspacing, size_template, grids = args.values()
    dimension = grids.shape[1]
    src = Template(self.source).render(
      gridspacing=gridspacing,
      dimension=dimension,
      maxdensity=1, #dummy
    )
    return str(src)
  def shape(self, args):
    coors, weights, origin, gridspacing, size_template, grids = args.values()
    assert np.dtype(coors.dtype) == np.dtype(weights.dtype) == np.float32
    assert len(coors.shape) == 3 and coors.shape[2] == 3, coors.shape
    assert len(weights.shape) == 1 and weights.shape[0] == coors.shape[1], (weights.shape, coors.shape)
    assert len(origin) == 3
    grids.dtype = np.float32 
    assert grids.shape[0] >= len(coors) #we are writing to a subset of the grids
    grids.shape = (len(coors),) + grids.shape[1:]
  def prepare(self, funcname, args):
    #__global__ void atomdensitymask1(
    # Coordinate *coor, float *weights, int size_coor, 
    # float *grids
    #)
    return GPUPrepare([np.intp, np.intp, np.int, np.intp])
  def define_grid(self, args):
    coors, weights, origin, gridspacing, size_template, grids = args.values()
    blocks_per_structure = int ( float(size_template) / self.THREADS_PER_BLOCK + 0.999 )  
    return GPUGrid((len(coors),blocks_per_structure ), (self.THREADS_PER_BLOCK,1,1))
  def kernel_args(self, args):
    return args["coors"], args["weights"], args["size_template"], args["grids"]


class atomdensitymask2(GPUAction):
  THREADS_PER_BLOCK = 256
  main = "atomdensitymask2"
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(output) == 1
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 5
    assert len(output) == 1
    d = OrderedDict(zip(("coors", "grids", "origin", "gridspacing", "maxdensity"),input_args))
    overlaps = output[0]
    d["size_template"] = d["coors"].shape[1]
    d["overlaps"] = overlaps
    return d
  def get_source(self, args):    
    coors, grids, origin, gridspacing, maxdensity, size_template, overlaps = args.values()
    dimension = grids.shape[1]
    src = Template(self.source).render(
      gridspacing=gridspacing,
      dimension=dimension,
      maxdensity=maxdensity,
    )
    return str(src)
  def shape(self, args):
    coors, grids, origin, gridspacing, maxdensity, size_template, overlaps = args.values()
    assert np.dtype(coors.dtype) == np.dtype(grids.dtype) == np.float32
    assert len(coors.shape) == 3 and coors.shape[2] == 3, coors.shape
    assert len(grids.shape) == 4 and grids.shape[0] == coors.shape[0], (grids.shape, coors.shape)
    assert grids.shape[1] == grids.shape[2] == grids.shape[3], grids.shape
    assert len(origin) == 3
    overlaps.dtype = np.float32 
    overlaps.shape = coors.shape[:2]
  def prepare(self, funcname, args):
    #__global__ void atomdensitymask2(
    # Coordinate *coor, float *overlaps, int size_coor, 
    # float *grids
    #)       
    return GPUPrepare([np.intp, np.intp, np.int, np.intp])
  def define_grid(self, args):
    coors, grids, origin, gridspacing, maxdensity, size_template, overlaps = args.values()
    blocks_per_structure = int ( float(size_template) / self.THREADS_PER_BLOCK + 0.999 )  
    return GPUGrid((len(coors),blocks_per_structure ), (self.THREADS_PER_BLOCK,1,1))
  def kernel_args(self, args):
    return args["coors"], args["overlaps"], args["size_template"], args["grids"]
