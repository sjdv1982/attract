"""
Calculate Gradient Vector Matching (GVM) score
usage: python gvm.py <SITUS map> <gradient threshold> <ATTRACT DAT file>  <ATTRACT-compatible PDB file>\
"""
import sys, os, copy
import multiprocessing

import numpy, scipy.signal, scipy.ndimage
import gridify
from math import sqrt

p = -1.0  #sign swaps, because we are *convoluting* with a kernel; an atom at X=10 means a negative X gradient at X=11!
gvm_kernel_z = numpy.array( [ [[-p, 0, p]] * 3] * 3)
gvm_kernel_x = gvm_kernel_z.swapaxes(2,0)
gvm_kernel_y = gvm_kernel_z.swapaxes(2,1)

def read_situs(situsfile):
  header = open(situsfile).readline()
  h = header.split()
  assert len(h) == 7, header
  voxelsize = float(h[0])
  origin = numpy.array([float(v) for v in h[1:4]])
  dimensions = tuple([int(v) for v in h[4:7]])
  data = numpy.loadtxt(situsfile,skiprows=2)
  data = data.reshape(-1)
  data = data.reshape(dimensions, order='F')
  return data, voxelsize, origin

def calc_gvm(coor):
  pdbmap = numpy.zeros(emdata.shape, order='F') 
  gridify.gridify(coor, pdbmap, gridspacing, origin)
  pdbgradx = scipy.ndimage.filters.prewitt(pdbmap, axis=0, mode="constant")[1:-1,1:-1,1:-1]
  pdbgrady = scipy.ndimage.filters.prewitt(pdbmap, axis=1, mode="constant")[1:-1,1:-1,1:-1]
  pdbgradz = scipy.ndimage.filters.prewitt(pdbmap, axis=2, mode="constant")[1:-1,1:-1,1:-1]
  sumx, sumxy, sumxx = 0.0,0.0,0.0
  for pdbgrad, emgradm,em_mask in zip((pdbgradx,pdbgrady,pdbgradz),(emgradxm,emgradym,emgradzm),(em_maskx,em_masky,em_maskz)):
    pdbgradm = pdbgrad[em_mask]
    sumx += numpy.sum(pdbgrad, axis=None)
    sumxx += numpy.sum(numpy.dot(pdbgradm, pdbgradm), axis=None)
    sumxy += numpy.sum(numpy.dot(pdbgradm, emgradm), axis=None)
  sxx = sumxx - sumx * sumx / corrcount;
  sxy = sumxy - sumx * sumy / corrcount;
  syy = sumyy - sumy * sumy / corrcount;
  return sxy/sqrt(sxx*syy)            
  
  
if __name__ == "__main__":    
  multiprocessing.freeze_support()
  import collectlibpy as collectlib

  ensfiles = []
  modefile = None
  imodefile = None

  anr = 0
  output = None
  while 1:
    anr += 1
        
    if anr > len(sys.argv)-1: break  
    arg = sys.argv[anr]
    
    if anr <= len(sys.argv)-3 and arg == "--ens":
      ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
      sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
      anr -= 3
      continue

    if anr <= len(sys.argv)-2 and arg == "--modes":
      modefile = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue
    
    if anr <= len(sys.argv)-2 and arg == "--imodes":
      imodefile = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue
    
    if arg.startswith("--"): raise Exception("Unknown option '%s'" % arg)
      
  emfile = sys.argv[1]
  assert os.path.exists(emfile), emfile
  threshold = float(sys.argv[2])
  datfile = sys.argv[3]
  assert os.path.exists(datfile), datfile
  pdbfile = sys.argv[4]
  assert os.path.exists(pdbfile), pdbfile
  
  emdata, gridspacing, origin = read_situs(emfile)

  emgradx = scipy.ndimage.filters.prewitt(emdata, axis=0, mode="constant")[1:-1,1:-1,1:-1]
  emgrady = scipy.ndimage.filters.prewitt(emdata, axis=1, mode="constant")[1:-1,1:-1,1:-1]
  emgradz = scipy.ndimage.filters.prewitt(emdata, axis=2, mode="constant")[1:-1,1:-1,1:-1]
  em_maskx = ((emgradx >= threshold) | (emgradx <= -threshold))
  em_masky = ((emgrady >= threshold) | (emgrady <= -threshold))
  em_maskz = ((emgradz >= threshold) | (emgradz <= -threshold))
  emgradxm = emgradx[em_maskx]
  emgradym = emgrady[em_masky]
  emgradzm = emgradz[em_maskz]
  
  sumy, sumyy = 0.0, 0.0
  corrcount = 0
  for emgradm, em_mask in zip((emgradxm,emgradym,emgradzm),(em_maskx,em_masky,em_maskz)):
    sumy += numpy.sum(emgradm, axis=None)
    sumyy += numpy.sum(emgradm * emgradm, axis=None)
    corrcount += numpy.count_nonzero(em_mask)
  
  initargs = [datfile, pdbfile]
  if modefile: initargs += ["--modes", modefile]
  if imodefile: initargs += ["--imodes", imodefile]
  for nr, ensfile in ensfiles:
    initargs += ["--ens", nr, ensfile]

  collectlib.collect_init(initargs)

  blocksize0 = 100
  poolsize = multiprocessing.cpu_count()
  blocksize = max(poolsize * int(blocksize0/poolsize), 1)
  pool = multiprocessing.Pool(poolsize)
  while 1:
    coors = []
    for nstruc in range(blocksize):
      eof = collectlib.collect_next()
      if eof: break    
      coor = collectlib.collect_coor(copy=True)          
      coors.append(coor)
    if eof and not nstruc: break
    #corr = [calc_gvm(coor) for coor in coors]
    corr = pool.map(calc_gvm, coors)
    for c in corr:
      print "%.6f" % c
    if eof: break  
