import sys, numpy, scipy.ndimage, scipy.signal
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"]+'/em')
from em_iolib import write_situs, write_mask
from math import *
from lib.pdb import get_coordinates
from scipy.spatial.distance import pdist

def make_mask(pdbfile,gridspacing,boxlim):
  coor = get_coordinates(pdbfile,ignore_weightless=True)
  com = numpy.mean(coor,axis=0)
  coor = coor-com
  coorx, coory, coorz = coor[:,0], coor[:,1], coor[:,2]
  subtractx = [numpy.abs(coorx-c)  for c in coorx]
  subtractx = numpy.reshape(subtractx,(-1))
  subtractx = subtractx[numpy.where(subtractx > 0)]
  subtracty = [numpy.abs(coory-c)  for c in coory]
  subtracty = numpy.reshape(subtracty,(-1))
  subtracty = subtracty[numpy.where(subtracty > 0)]
  subtractz = [numpy.abs(coorz-c)  for c in coorz]
  subtractz = numpy.reshape(subtractz,(-1))
  subtractz = subtractz[numpy.where(subtractz > 0)]
  pdbgridspacing = [numpy.min(subtractx),numpy.min(subtracty),numpy.min(subtractz)]
  assert max(pdbgridspacing) <= 2.0*gridspacing
  #assert numpy.amax(coor) < boxlim, numpy.amin(coor) > -1.0*boxlim
  dimensions = boxlim/gridspacing*2.0
  #x increment change fastest in situs format, z increments slowest
  grid = [0 for x in range(int(dimensions**3.0))]
  dimensions = (dimensions,dimensions,dimensions)
  grid = numpy.reshape(grid,dimensions,order='F')
  origin = -boxlim
  for xx, yy, zz in coor:
    resetx = (xx - origin)/gridspacing
    resety = (yy - origin)/gridspacing
    resetz = (zz - origin)/gridspacing
    resetx = int(resetx)
    resety = int(resety)
    resetz = int(resetz)
    if resetx >= dimensions[0] or resety >= dimensions[0] or resetz >= dimensions[0] or resetx < 0 or resety < 0 or resetz < 0: continue 
    #print xx, yy, zz
    #print "sorted in", resetx, resety, resetz
    #print "voxel coordinates", resetx*gridx+originx,resety*gridy+originy,resetz*gridz+originz
    #print "gridspacing", gridx, gridy,gridz
    #print "origin" , origin
    #sys.exit()
    grid[resetx][resety][resetz] = 1.0

  return grid

try:
  pdbfile = sys.argv[1] #dummy residue low resolution model from DAMMFILT, centered at origin
  gridspacing = float(sys.argv[2])
  boxlim = float(sys.argv[3]) #construct a mask from -boxlim to boxlim Angstroms, in all directions
  assert boxlim/gridspacing <= 20 #limit for the GPU 40x40x40 grid
  maskfile = sys.argv[4]
  outputsitusfile = None
  if len(sys.argv) == 6:
    outputsitusfile = sys.argv[5] 
except:
  import traceback
  traceback.print_exc()
  print
  print("Syntax: pdb2mask <dummy residue pdb file> <mask voxelsize> <box limit in A> <output mask file in .mask format> [output mask file in SITUS format]")
  sys.exit(1)
  
extrusion = 5 #2 on each side! changed from cryoEM, yields larger margin 
origin = (0,0,0)    
mask = make_mask(pdbfile,gridspacing,boxlim)
print mask.shape
extrusionkernel = numpy.array([1.0]*extrusion**3).reshape((extrusion,extrusion,extrusion))
mask = scipy.signal.convolve(mask, extrusionkernel, "same")
mask[mask>=0.02]=1
boolmask = numpy.greater(mask, 0)
print boolmask[:10]
maskorigin = (-boxlim,-boxlim,-boxlim)
pad = []
maskpadorigin = []
print mask.shape
#for ori,d in zip(maskorigin, mask.shape):
  #lpad = int(ceil(ori+boxlim)/gridspacing)
  #end = ori + gridspacing * d
  #rpad = int(ceil(boxlim-end)/gridspacing)
  #pad.append((lpad, lpad+d, lpad+rpad+d))
  #print "padding", lpad, lpad+d,lpad+rpad+d
  #maskpadorigin.append(ori-lpad*gridspacing)
maskpadorigin=numpy.array(maskpadorigin)  
#padmask = numpy.zeros([v[2] for v in pad], dtype=float,order='F')
#print padmask.shape
#padboolmask = numpy.zeros([v[2] for v in pad], dtype=bool,order='F')
#print padboolmask.shape
#print "pad",len(pad), pad
#padboolmask[pad[0][0]:pad[0][1],pad[1][0]:pad[1][1],pad[2][0]:pad[2][1]] = boolmask
#padmask[pad[0][0]:pad[0][1],pad[1][0]:pad[1][1],pad[2][0]:pad[2][1]] = mask

write_mask(maskfile, boolmask, gridspacing, maskorigin)

if outputsitusfile:
  write_situs(outputsitusfile, mask, gridspacing, maskorigin )


