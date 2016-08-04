from topview import GPUAction, TVArray
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os
from jinja2 import Template

sourcefile = os.path.split(__file__)[0] + os.sep + "/gvm.cu"

class _gvm_base(GPUAction):
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 2
    assert len(output) == 3
    d = OrderedDict(zip(("maps", "refe_map"),input_args))
    sumx, sumxx, sumxy = output
    d["sumx"] = sumx
    d["sumxx"] = sumxx
    d["sumxy"] = sumxy
    return d
  def get_source(self, args):    
    maps, refe_map, sumx, sumxx, sumxy = args.values()
    src = Template(self.source).render(
      dimensions=maps.shape[1:],
      dimensions2=refe_map.shape,
    )
    return str(src)
  def shape(self, args):
    maps, refe_map, sumx, sumxx, sumxy = args.values()
    assert np.dtype(maps.dtype) == np.dtype(refe_map.dtype) == np.float32
    assert len(refe_map.shape) == 3, refe_map.shape
    assert refe_map.shape[0] == maps.shape[1]-2, (refe_map.shape, maps.shape)
    assert refe_map.shape[1] == maps.shape[2]-2, (refe_map.shape, maps.shape)
    assert refe_map.shape[2] == maps.shape[3]-2, (refe_map.shape, maps.shape)
    sumx.dtype = np.float32 
    sumx.shape = (len(maps),)
    sumxx.dtype = np.float32 
    sumxx.shape = (len(maps),)
    sumxy.dtype = np.float32 
    sumxy.shape = (len(maps),)
  def prepare(self, funcname, args):
    #__global__ void gvm_x(Map *maps, RefeMap *refe_map0, float *sumx, float *sumxx, float *sumxy) {
    return GPUPrepare([np.intp, np.intp, np.intp, np.intp, np.intp])
  def define_grid(self, args):
    maps, refe_map, sumx, sumxx, sumxy = args.values()
    d = maps.shape[1:]
    blocky = int((d[0] - 2)/30.0 + 0.99999)
    blockz = d[self.blockdim] - 2
    blockyz = blocky * blockz
    return GPUGrid((len(maps), blockyz), (32,1,1))
  def kernel_args(self, args):
    return args["maps"], args["refe_map"], args["sumx"], args["sumxx"], args["sumxy"]

class gvm_x(_gvm_base):
  main = "gvm_x"
  blockdim = 2

class gvm_y(_gvm_base):
  main = "gvm_y"
  blockdim = 1

class gvm_z(_gvm_base):
  main = "gvm_z"
  blockdim = 2
  