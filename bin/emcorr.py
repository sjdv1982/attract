"""
Calculate EM cross-correlation score
usage: python emcorr.py <SITUS map> <resolution> <threshold> <ATTRACT DAT file>  <ATTRACT-compatible PDB file>\
"""
import sys, os, copy
import multiprocessing

import numpy, scipy.signal, scipy.ndimage
import gridify
from math import sqrt

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

def calc_emcorr(coor):
  pdbmap = numpy.zeros(emdata.shape, order='F') 
  gridify.gridify(coor, pdbmap, gridspacing, origin)
  pdbmap_gaussian = scipy.ndimage.filters.gaussian_filter(pdbmap, sigma, mode="constant")
  pdb_mask = (pdbmap_gaussian > threshold)
  corrcount = numpy.count_nonzero(pdb_mask)
  pdbmap_gaussianm = pdbmap_gaussian[pdb_mask]
  emdata_gaussianm = emdata_gaussian[pdb_mask]
  sumy = numpy.sum(emdata_gaussianm, axis=None)
  sumyy = numpy.sum(emdata_gaussianm * emdata_gaussianm, axis=None)  
  sumx = numpy.sum(pdbmap_gaussianm, axis=None)
  sumxx = numpy.sum(numpy.dot(pdbmap_gaussianm, pdbmap_gaussianm), axis=None)
  sumxy = numpy.sum(numpy.dot(pdbmap_gaussianm, emdata_gaussianm), axis=None)
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
  resolution = float(sys.argv[2])
  threshold = float(sys.argv[3])
  datfile = sys.argv[4]
  assert os.path.exists(datfile), datfile
  pdbfile = sys.argv[5]
  assert os.path.exists(pdbfile), pdbfile
  
  emdata, gridspacing, origin = read_situs(emfile)
  sig1 = resolution/2.0;
  varmap = sig1*sig1/(gridspacing*gridspacing); 
  sigma = sqrt(varmap/3.0); 
  emdata_gaussian = scipy.ndimage.filters.gaussian_filter(emdata, sigma, mode="constant")
    
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
    #corr = [calc_emcorr(coor) for coor in coors]
    corr = pool.map(calc_emcorr, coors)
    for c in corr:
      print "%.6f" % c
    if eof: break  
