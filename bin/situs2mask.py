import sys, numpy, scipy.ndimage, scipy.signal, struct
from math import *

try:
  situsfile = sys.argv[1]
  threshold = float(sys.argv[2])
  maskvoxelsize = float(sys.argv[3])
  maskfile = sys.argv[4]
  outputsitusfile = None
  if len(sys.argv) == 6:
    outputsitusfile = sys.argv[5] 
except:
  import traceback
  traceback.print_exc()
  print
  print("Syntax: situs2mask <SITUS map file> <map threshold> <mask voxelsize> <output mask file in .mask format> [output mask file in SITUS format]")
  sys.exit(1)
  
  
def read_situs(situsfile):
  header = open(situsfile).readline()
  h = header.split()
  assert len(h) == 7, header
  voxelsize = float(h[0])
  origin = numpy.array([float(v) for v in h[1:4]])
  dimensions = tuple([int(v) for v in h[4:7]])
  data = numpy.genfromtxt(situsfile,skip_header=2,skip_footer=1)
  nvox = dimensions[0]*dimensions[1]*dimensions[2]
  for lastline in open(situsfile).xreadlines():
    pass
  lastdata = numpy.array([float(v) for v in lastline.split()])
  data = data.reshape(-1)
  data = numpy.append(data, lastdata)
  data = data.reshape(dimensions, order='F')
  return data, voxelsize, origin


def write_situs(situsfile, data, voxelsize, origin):
  f = open(situsfile, "w")
  f.write("%.6f " % voxelsize)
  f.write("%.6f %.6f %.6f " % tuple(origin))
  f.write("%d %d %d\n\n" % data.shape)
  d = data.reshape(-1, order='F')
  count = 0
  for v in d:
    f.write("%11.6f " % v)
    count += 1
    if not count % 10: f.write("\n")
    #else: f.write(" ")
  f.close()

def write_mask(maskfile, data, voxelsize, origin):
  f = open(maskfile, "wb") 
  f.write("ATTRACTMASK")  
  voxelsize = struct.pack("f", voxelsize)
  f.write(voxelsize)  
  origin = struct.pack("fff", *tuple(origin))
  f.write(origin)  
  shape = struct.pack("iii", *data.shape)
  f.write(shape)
  f.write(data.data)
  
  
data, gridspacing, origin = read_situs(situsfile)
data[data<threshold]=-1
data[data>=threshold]=1
assert gridspacing < maskvoxelsize
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
extrusion = numpy.array([1/27.0]*27).reshape((3,3,3))
mask = scipy.signal.convolve(mask, extrusion, "same")
mask[mask>=0.02]=1
boolmask = numpy.greater(mask, 0)
originshift = (ratio-1)*gridspacing
maskorigin = origin + originshift
write_mask(maskfile, boolmask, maskvoxelsize, maskorigin )

if outputsitusfile:
  write_situs(outputsitusfile, mask, maskvoxelsize, maskorigin )


