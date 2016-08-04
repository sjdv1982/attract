from topview import GPUAction, TVArray
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os
from jinja2 import Template

sourcefile = os.path.split(__file__)[0] + os.sep + "/calc_radist.cu"

class calc_radist(GPUAction):
  THREADS_PER_BLOCK = 1024
  main = "calc_radist"
  source = open(sourcefile).read()
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert accumulate
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 8
    assert len(output) == 1
    coors1, vacweights1, dummyweights1, coors2, vacweights2,dummyweights2, binsize, nbins = input_args

    assert len(coors1.shape) == len(coors2.shape) == 3
    assert coors1.shape[2] == coors2.shape[2] == 3
    assert len(vacweights1.shape) == len(vacweights2.shape) == 1
    assert len(dummyweights1.shape) == len(dummyweights2.shape) == 1
    assert coors1.shape[0] == coors2.shape[0]
    nstruc = coors1.shape[0]
    assert coors1.shape[1] == len(vacweights1)
    assert coors2.shape[1] == len(vacweights2)
    assert coors1.shape[1] == len(dummyweights1)
    assert coors2.shape[1] == len(dummyweights2)
    coorlen = coors1.shape[1], coors2.shape[1]
    binsize = float(binsize)
    nbins = int(nbins)
    
    radist = output[0]
    assert radist.shape is not None and len(radist.shape) == 3 #radist must be pre-shaped
    assert radist.shape[0] == nstruc
    assert radist.shape[1] == 3
    assert radist.shape[2] == nbins
    
    for v in "coors1", "vacweights1", "dummyweights1", "coors2", "vacweights2", "dummyweights2", "binsize", "nbins", "radist", "coorlen", "nstruc":
      d[v] = locals()[v]
      
    return d
  def get_source(self, args):
    coors1, vacweights1, dummyweights1, coors2, vacweights2, dummyweights2, binsize, nbins, radist, coorlen, nstruc = args.values()
    src = Template(self.source).render(
      coorlen = coorlen,
      nstruc=nstruc,
      nbins=nbins,
      binsize=binsize,
      THREADS_PER_BLOCK=self.THREADS_PER_BLOCK,
    )
    return str(src)
  def shape(self, args):
    coors1, vacweights1, dummyweights1, coors2, vacweights2, dummyweights2, binsize, nbins, radist, coorlen, nstruc = args.values()
    assert np.dtype(coors1.dtype) == np.float32
    assert np.dtype(coors2.dtype) == np.float32
    assert np.dtype(vacweights1.dtype) == np.float32
    assert np.dtype(vacweights2.dtype) == np.float32
    assert np.dtype(dummyweights1.dtype) == np.float32
    assert np.dtype(dummyweights2.dtype) == np.float32
    radist.dtype = np.float32
  def prepare(self, funcname, args):
    """
    __global__ void calc_radist(
    Coordinates1 *coors1, float *weights1,
    Coordinates2 *coors2, float *weights2,
    Radists *radist
    ) {  
    """  
    return GPUPrepare([np.intp, np.intp, np.intp, np.intp, np.intp, np.intp, np.intp])
  def define_grid(self, args):
    coors1, vacweights1, dummyweights1, coors2, vacweights2, dummyweights2, binsize, nbins, radist, coorlen, nstruc = args.values()
    blocks_per_structure = int ( float(coorlen[0]) / self.THREADS_PER_BLOCK + 0.999 )  
    return GPUGrid((nstruc,blocks_per_structure ), (self.THREADS_PER_BLOCK,1,1))
  
  def kernel_args(self, args):
    return args["coors1"], args["vacweights1"], args["dummyweights1"], args["coors2"], args["vacweights2"], args["dummyweights2"], args["radist"]
