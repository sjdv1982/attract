from topview import GPUAction
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os

sourcefile = os.path.split(__file__)[0] + os.sep + "/euler2rotmat.cu"

class euler2rotmat(GPUAction):
  THREADS_PER_BLOCK = 256
  main = "euler2rotmat"
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(output) == 1
    d = OrderedDict()
    assert len(input_args) + len(input_kwargs) == 1
    if len(input_args) >= 1:
      d["eulers"] = input_args[0]
    else:
      d["eulers"] = input_kwargs["eulers"]          
    d["rotmats"] = output[0]
    d["arraylength"] = len(d["eulers"])
    return d
  def shape(self, args):
    eulers, rotmats, arraylength = args.values()
    assert eulers.dtype == np.float32    
    assert len(eulers.shape) == 2 and eulers.shape[1] == 3, eulers.shape
    rotmats.dtype = np.float32
    rotmats.shape = (eulers.shape[0], 9)
  def prepare(self, funcname, args):
    return GPUPrepare([np.intp, np.intp, np.int])
  def define_grid(self, args):
    eulers = args["eulers"]
    blocks = int ( float(len(eulers)) / self.THREADS_PER_BLOCK + 0.999 )
    return GPUGrid((blocks,1), (self.THREADS_PER_BLOCK,1,1))
