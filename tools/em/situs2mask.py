import sys, numpy, scipy.ndimage, scipy.signal
from em_iolib import read_situs, write_situs, write_mask
from math import *

try:
  situsfile = sys.argv[1]
  threshold = float(sys.argv[2])
  maskvoxelsize = float(sys.argv[3])
  boxlim = float(sys.argv[4]) #construct a mask from -boxlim to boxlim Angstroms, in all directions
  maskfile = sys.argv[5]
  outputsitusfile = None
  if len(sys.argv) == 7:
    outputsitusfile = sys.argv[6] 
except:
  import traceback
  traceback.print_exc()
  print
  print("Syntax: situs2mask <SITUS map file> <map threshold> <mask voxelsize> <box limit in A> <output mask file in .mask format> [output mask file in SITUS format]")
  sys.exit(1)
  
extrusion = 3 #1 on each side
    
data, gridspacing, origin = read_situs(situsfile)
data[data<threshold]=-1
data[data>=threshold]=1
assert gridspacing <= maskvoxelsize
ratio = maskvoxelsize/gridspacing
ceilratio = int(ceil(ratio))
zoomratio = ceilratio/ratio
zoomdata = scipy.ndimage.interpolation.zoom(data, zoomratio)
padzoomdatasize = [ceilratio*int(ceil(v/float(ceilratio))) for v in zoomdata.shape]
padzoomdata = numpy.zeros(padzoomdatasize,dtype=zoomdata.dtype)
padzoomdata[:zoomdata.shape[0], :zoomdata.shape[1], :zoomdata.shape[2]] = zoomdata

masksize = [v/ceilratio for v in padzoomdata.shape]
mask = padzoomdata.reshape([masksize[0], ceilratio, masksize[1], ceilratio, masksize[2], ceilratio]).mean(5).mean(3).mean(1)
mask[mask>=0]=1
mask[mask<0]=0
extrusionkernel = numpy.array([1.0]*extrusion**3).reshape((extrusion,extrusion,extrusion))
mask = scipy.signal.convolve(mask, extrusionkernel, "same")
mask[mask>=0.02]=1
boolmask = numpy.greater(mask, 0)
originshift = (ratio-1)*gridspacing
maskorigin = origin + originshift
pad = []
maskpadorigin = []
for ori,d in zip(maskorigin, mask.shape):
  lpad = int(ceil(ori+boxlim)/maskvoxelsize)
  end = ori + maskvoxelsize * d
  rpad = int(ceil(boxlim-end)/maskvoxelsize)
  pad.append((lpad, lpad+d, lpad+rpad+d))
  maskpadorigin.append(ori-lpad*maskvoxelsize)
maskpadorigin=numpy.array(maskpadorigin)  
padmask = numpy.zeros([v[2] for v in pad], dtype=float,order='F')
padboolmask = numpy.zeros([v[2] for v in pad], dtype=bool,order='F')
padboolmask[pad[0][0]:pad[0][1],pad[1][0]:pad[1][1],pad[2][0]:pad[2][1]] = boolmask
padmask[pad[0][0]:pad[0][1],pad[1][0]:pad[1][1],pad[2][0]:pad[2][1]] = mask

write_mask(maskfile, padboolmask, maskvoxelsize, maskpadorigin)

if outputsitusfile:
  write_situs(outputsitusfile, padmask, maskvoxelsize, maskpadorigin )


