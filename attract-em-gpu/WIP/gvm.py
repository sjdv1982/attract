from lib.read_struc import read_struc, depivotize
from lib.read_situs import read_situs
from lib.pdb import get_coordinates_multi_attract, get_weights
from lib.euler2rotmat import euler2rotmat
from lib.collect import collect
from lib.gridify import gridify
from lib.core import copy, transpose, fill, nptask
from lib.atomdensitymask import atomdensitymask1, atomdensitymask2
from lib.gvm import gvm_x, gvm_y, gvm_z
import numpy as np
import scipy.signal, scipy.ndimage

from topview import TVArray, TVContext, TVArrayList, TVHash

grid_chunksize = 1000

import sys
emfile = sys.argv[1]
threshold = float(sys.argv[2])
datfile = sys.argv[3]
pdb = sys.argv[4]


pivots, eulers, positions0 = read_struc(datfile, None) #TODO auto-pivots
positions = depivotize(pivots, eulers, positions0)

emdata, gridspacing, origin = read_situs(emfile)
#emdata = emdata[10:17,10:17,10:17]

templates = get_coordinates_multi_attract(pdb, ignore_weightless=False)
#weights = [get_weights(pdb, ignore_weightless=False) for pdb in pdbfiles]
#weights = [1.0 for w in weights] ### original gvm.py does not weight atoms

emdata = np.ascontiguousarray(emdata, dtype=np.float32)
emgrads = []
sumy = 0
sumyy = 0
corrcount = 0
for n in range(3):
  emgrad = scipy.ndimage.filters.prewitt(emdata, axis=n, mode="constant")[1:-1,1:-1,1:-1]
  emgrad = np.ascontiguousarray(emgrad, dtype=np.float32)
  emgrad[(emgrad > -threshold) & (emgrad < threshold)] = 0
  dsumy, dsumyy = emgrad.sum(), (emgrad*emgrad).sum()
  sumy += dsumy
  sumyy += dsumyy
  corrcount += np.count_nonzero(emgrad)
  emgrads.append(emgrad)
syy = sumyy - sumy * sumy / corrcount;
emgrads *= 3 ###

def gvm(refe, eulers, positions, templates, gridshape, origin, gridspacing, chunksize):
  chunklen = min(chunksize, len(eulers))
  
  nbodies = len(templates)
  assert eulers.shape[1] == positions.shape[1] == nbodies
  nstruc = eulers.shape[0]
  assert eulers.shape[0] == positions.shape[0]
  
  e = TVArray("e", eulers)  
  p = TVArray("p", positions)
  eT = TVArray("eT")  
  pT = TVArray("pT")
  g_e = TVArray("g_e", gpu=True)  
  g_p = TVArray("g_p", gpu=True)
  g_rotmats = TVArray("g_rotmats", gpu=True)
  g_coors = []
  g_templates = []
  i_templates = TVArrayList("i_templates", templates)
  i_refe = TVArrayList("i_refe", refe)
  for b in range(nbodies):
    g_coors.append(TVArray("g_coors[%d]" % b, gpu=True))
    g_templates.append(TVArray("g_templates[%d]" % b, gpu=True))
    copy(i_templates[b]) > g_templates[b]
  transpose(e, (1,0,2)) > eT
  transpose(p, (1,0,2)) > pT
  copy(eT) > g_e
  copy(pT) > g_p
  g_ee = g_e.rsplit()
  g_pp = g_p.rsplit()

  for b in range(nbodies):
    euler2rotmat(g_ee[b]) > g_rotmats      
    collect(g_pp[b], g_rotmats, g_templates[b]) > g_coors[b]
  
  chunk_coors = []
  for b in range(nbodies):
    chunk_coors.append(g_coors[b].rchunks(chunklen))
  
  maps = TVArray("maps", gpu=True, dtype="float32", shape = (chunklen,) + tuple(gridshape))
  g_refe = []
  for n in range(3):
    assert gridshape[0] == refe[n].shape[0] + 2
    assert gridshape[1] == refe[n].shape[1] + 2
    assert gridshape[2] == refe[n].shape[2] + 2
    g_refe.append(TVArray("g_refe[%d]" % n, gpu=True))
    copy(i_refe[n]) > g_refe[n]
  g_sumx, g_sumxx, g_sumxy = [],[],[]
  sumx, sumxx, sumxy = [],[],[]
  chunk_sumx, chunk_sumxx, chunk_sumxy = [],[],[]
  for n in range(3):
    g_sumx.append(TVArray("g_sumx[%d]" % n, gpu=True, shape = (chunklen,), dtype="float32" ))
    g_sumxx.append(TVArray("g_sumxx[%d]" % n, gpu=True, shape = (chunklen,), dtype="float32" ))
    g_sumxy.append(TVArray("g_sumxy[%d]" % n, gpu=True, shape = (chunklen,), dtype="float32" ))
    sumx.append(TVArray("sumx[%d]" % n, shape = (len(eulers),), dtype="float32" ))
    sumxx.append(TVArray("sumxx[%d]" % n, shape = (len(eulers),), dtype="float32" ))
    sumxy.append(TVArray("sumxy[%d]" % n, shape = (len(eulers),), dtype="float32" ))
    chunk_sumx.append(sumx[n].wchunks(chunksize))
    chunk_sumxx.append(sumxx[n].wchunks(chunksize))
    chunk_sumxy.append(sumxy[n].wchunks(chunksize))
  
  for i in range(len(chunk_sumx[0])):  
    fill(0) > maps
    for b in range(nbodies):
      gridify(chunk_coors[b][i], origin, gridspacing) >> maps
    gvm_x(maps, g_refe[0]) > (g_sumx[0], g_sumxx[0], g_sumxy[0])
    gvm_y(maps, g_refe[1]) > (g_sumx[1], g_sumxx[1], g_sumxy[1])                       
    gvm_z(maps, g_refe[2]) > (g_sumx[2], g_sumxx[2], g_sumxy[2])
    for n in range(3):
      copy(g_sumx[n]) > chunk_sumx[n][i]
      copy(g_sumxx[n]) > chunk_sumxx[n][i]
      copy(g_sumxy[n]) > chunk_sumxy[n][i]
  for n in range(3):
    sumx[n].join()
    sumxx[n].join()
    sumxy[n].join()
  return sumx[0], sumx[1], sumx[2], sumxx[0], sumxx[1], sumxx[2], sumxy[0], sumxy[1], sumxy[2]

ret = gvm(emgrads, eulers, positions, templates, tuple(emdata.shape), tuple(origin), gridspacing, grid_chunksize)
sumx = list(ret[0:3])
sumxx = list(ret[3:6])
sumxy = list(ret[6:9])
for n in range(3):
  sumx[n] = sumx[n].get_data()
  sumxx[n] = sumxx[n].get_data()
  sumxy[n] = sumxy[n].get_data()

"""
from  topview import TVContext
gvm  = TVContext(gvm)
ret = gvm(emgrads, eulers, positions, templates, tuple(emdata.shape), tuple(origin), gridspacing, grid_chunksize)
sumx = list(ret[0:3])
sumxx = list(ret[3:6])
sumxy = list(ret[6:9])
"""

for s in range(len(sumx[0])):
  csumxx, csumx, csumxy = 0,0,0
  for n in range(3):
    csumx += sumx[n][s]
    csumxx += sumxx[n][s]
    csumxy += sumxy[n][s]
  sxx = csumxx - csumx * csumx / corrcount;
  sxy = csumxy - csumx * sumy / corrcount;  
  variance = sxx*syy
  corr = 0.0
  if variance > 0.0000001: 
    corr = sxy/np.sqrt(variance)            
  print "%.5f" % corr
    
  