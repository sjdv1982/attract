from topview import GPUAction
from topview.TVAction import GPUPrepare, GPUGrid
from collections import OrderedDict
import numpy as np
import os

sourcefile = os.path.split(__file__)[0] + os.sep + "/collect.cu"

class collect(GPUAction):
  THREADS_PER_BLOCK = 256
  main = "collect"
  source = open(sourcefile).read()
  
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(output) == 1
    d = OrderedDict()
    assert len(input_kwargs) == 0
    assert len(input_args) == 3
    assert len(output) == 1
    positions, rotmats, template = input_args
    d = OrderedDict(zip(("rotmats", "translations", "template"),(rotmats, positions, template)))
    d["size_template"] = len(template)
    d["output"] = output[0]
    return d
  def shape(self, args):
    rotmats, translations, template, size_template, output = args.values()
    assert rotmats.dtype == translations.dtype == template.dtype == np.float32
    assert len(rotmats.shape) == 2 and rotmats.shape[1] == 9, rotmats.shape
    assert len(translations.shape) == 2 and translations.shape[1] == 3, translations.shape
    assert translations.shape[0] == rotmats.shape[0]
    assert len(template.shape) == 2 and template.shape[1] == 3, translations.shape
    output.dtype = np.float32
    oshape = translations.shape[:1] + template.shape
    output.shape = oshape
    assert output.shape == oshape
  def prepare(self, funcname, args):
     #__global__ void collect(RotMat *rotmats, Coordinate *translations, Coordinate *tmplate, int size_template, Coordinate *coordinates)
    return GPUPrepare([np.intp, np.intp, np.intp, np.int, np.intp])
  def define_grid(self, args):    
    size_template = args["size_template"]
    translations = args["translations"]
    blocks_per_structure = int ( float(size_template) / self.THREADS_PER_BLOCK + 0.999 )  
    return GPUGrid((len(translations),blocks_per_structure ), (self.THREADS_PER_BLOCK,1,1))
    