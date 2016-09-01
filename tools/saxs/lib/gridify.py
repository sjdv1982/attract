from topview import GPUAction, TVArray
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os
from jinja2 import Template

sourcefile = os.path.split(__file__)[0] + os.sep + "/gridify.cu"

class gridify(GPUAction):
  THREADS_PER_BLOCK = 256
  main = "gridify"
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert accumulate
    assert len(output) == 1
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 3
    assert len(output) == 1
    d = OrderedDict(zip(("coors", "origin", "gridspacing"),input_args))
    grids = output[0]
    assert grids.shape is not None and len(grids.shape) == 4 #grids must be pre-shaped
    d["size_template"] = d["coors"].shape[1]
    d["dim"] = grids.shape[1:]
    d["grids"] = grids
    return d
  def get_source(self, args):
    coors, origin, gridspacing, size_template, dim, grids = args.values()
    src = Template(self.source).render(
      origin=origin,
      gridspacing=gridspacing,
      dim=dim,
    )
    return str(src)
  def shape(self, args):
    coors, origin, gridspacing, size_template, dim, grids = args.values()
    assert np.dtype(coors.dtype) == np.float32
    assert len(coors.shape) == 3 and coors.shape[2] == 3, coors.shape
    assert len(origin) == 3
    grids.dtype = np.float32 
    assert grids.shape[0] >= len(coors) #we are writing to a subset of the grids
    grids.shape = (len(coors),) + grids.shape[1:]
  def prepare(self, funcname, args):
    #__global__ void gridify(
    #  Coordinate *coor, int size_coor, 
    #  float *grids
    #)     
    return GPUPrepare([np.intp, np.int, np.intp])
  def define_grid(self, args):
    coors, origin, gridspacing, size_template, dim, grids = args.values()
    blocks_per_structure = int ( float(size_template) / self.THREADS_PER_BLOCK + 0.999 )  
    return GPUGrid((len(coors),blocks_per_structure ), (self.THREADS_PER_BLOCK,1,1))

  """
  def convert(self, args):
    for k,v in args.items():
      if isinstance(v, TVArray):
        print k, v.dtype, v.shape
    return GPUAction.convert(self, args)
  """
  
  def kernel_args(self, args):
    return args["coors"], args["size_template"], args["grids"]
