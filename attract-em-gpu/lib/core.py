from topview import TVAction, TVArray
from collections import OrderedDict
import numpy as np

class copy(TVAction):
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(input_kwargs) == 0
    assert len(input_args) == 1
    assert len(output) == 1
    i, o = input_args + output
    i_gpu = False
    if i.gpu: i_gpu = True
    o_gpu = False
    if o.gpu: o_gpu = True    
    return OrderedDict(zip(("i", "o", "i_gpu", "o_gpu"),(i, o, i_gpu, o_gpu)))
  def shape(self, args):
    i, o, i_gpu, o_gpu = args.values()
    assert isinstance(i, TVArray), i
    assert isinstance(o, TVArray), o
    o.dtype = i.dtype
    o.shape = i.shape
  def run(self, args):
     i, o, i_gpu, o_gpu = args.values()
     if not i_gpu and not o_gpu:
       o[:] = i
     elif i_gpu and not o_gpu:
       i.get(o)
     elif o_gpu and not i_gpu:
       o.set(i)
     elif o_gpu and i_gpu:  
       #i.get(o)
       import pycuda
       assert i.nbytes == o.nbytes
       pycuda.driver.memcpy_dtod(o.gpudata, i.gpudata, i.nbytes)

class transpose(TVAction):
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(input_kwargs) == 0
    assert len(input_args) == 2
    assert len(output) == 1
    inp, axes, outp = input_args + output
    return OrderedDict(zip(("inp", "axes", "outp"),(inp, axes, outp)))
  def shape(self, args):
    inp, axes, outp = args.values()
    assert not inp.gpu
    assert not outp.gpu
    assert isinstance(inp, TVArray), inp
    assert isinstance(outp, TVArray), outp
    assert len(axes) == len(inp.shape)
    outp.dtype = inp.dtype
    oshape = []
    for n in range(len(axes)):
      oshape.append(inp.shape[axes[n]])      
    outp.shape = tuple(oshape)
  def run(self, args):
    inp, axes, outp, = args.values()
    outp[:] = np.ascontiguousarray(inp.transpose(axes))

class fill(TVAction):
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert not accumulate
    assert len(input_kwargs) == 0
    assert len(input_args) == 1
    assert len(output) == 1
    value, a = input_args + output
    return OrderedDict(zip(("a", "value"),(a,value)))
  def shape(self, args):
    a, value = args.values()
    assert a.dtype is not None
    assert a.shape is not None
  def run(self, args):
    a, value = args.values()    
    a.fill(value)


class nptask(TVAction):
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    assert len(input_kwargs) == 0
    assert len(input_args) == 2
    assert len(output) == 1
    inp, expr, output = input_args + output
    return OrderedDict(zip(("inp", "expr", "output", "accumulate"),(inp, expr, output, accumulate)))
  def shape(self, args):
    inp, expr, output, accumulate = args.values()
    assert output.dtype is not None #we cannot predict the dtype of the expression
    assert output.shape is not None #we cannot predict the shape of the expression
  def run(self, args):
    inp, expr, output, accumulate = args.values()
    v = eval("inp" + expr)
    if accumulate:
      output[:] += v
    else:
      output[:] = v
    
