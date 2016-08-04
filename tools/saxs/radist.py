from lib.read_struc import read_struc, depivotize
from lib.pdb import get_coordinates, get_saxs_factors
from lib.euler2rotmat import euler2rotmat
from lib.collect import collect
from lib.core import copy, transpose, fill
from lib.calc_radist import calc_radist
import numpy as np

from topview import TVArray, TVContext, TVArrayList

import sys
datfile = sys.argv[1]
pdbs = sys.argv[2:]

pivots, eulers, positions0 = read_struc(datfile, None) #TODO auto-pivots
positions = depivotize(pivots, eulers, positions0)
templates = [get_coordinates(pdb, ignore_weightless=False) for pdb in pdbs]
saxs_factors = [get_saxs_factors(pdb) for pdb in pdbs]

CHUNKSIZE = 2000

def func_collect(eulers, positions, templates):
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
  o_coors = []
  g_templates = []
  i_templates = TVArrayList("i_templates", templates)
  for b in range(nbodies):
    g_coors.append(TVArray("g_coors[%d]" % b, gpu=True))
    o_coors.append(TVArray("o_coors[%d]" % b))
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
    copy(g_coors[b]) > o_coors[b]
  
  return tuple([o for o in o_coors])
                   
def func_radist(eulers, positions, templates, saxs_factors, binsize, nbins):
  
  nbodies = len(templates)
  assert eulers.shape[1] == positions.shape[1] == nbodies
  nstruc = eulers.shape[0]
  assert eulers.shape[0] == positions.shape[0]

  o_radist = TVArray("radist", shape = (nstruc, nbins), dtype="float32") 
  
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
  g_saxs_factors = []
  i_saxs_factors = TVArrayList("i_saxs_factors", saxs_factors)  
  for b in range(nbodies):
    g_coors.append(TVArray("g_coors[%d]" % b, gpu=True))
    g_templates.append(TVArray("g_templates[%d]" % b, gpu=True))
    copy(i_templates[b]) > g_templates[b]
    g_saxs_factors.append(TVArray("g_saxs_factors[%d]" % b, gpu=True))
    copy(i_saxs_factors[b]) > g_saxs_factors[b]    
  transpose(e, (1,0,2)) > eT
  transpose(p, (1,0,2)) > pT
  copy(eT) > g_e
  copy(pT) > g_p
  g_ee = g_e.rsplit()
  g_pp = g_p.rsplit()

  for b in range(nbodies):
    euler2rotmat(g_ee[b]) > g_rotmats      
    collect(g_pp[b], g_rotmats, g_templates[b]) > g_coors[b]
  
  
  coors_chunks = [g.rchunks(CHUNKSIZE) for g in g_coors]
  o_radist_chunks = o_radist.wchunks(CHUNKSIZE)  
  g_radist = TVArray("g_radist", gpu = True, shape = (CHUNKSIZE, nbins), dtype="float32")  
  for i in range(len(o_radist_chunks)):
    g_radist.shape = (coors_chunks[b][i].shape[0], nbins)
    fill(0) > g_radist
    for b in range(nbodies):
      for bb in range(b+1, nbodies):
        calc_radist(coors_chunks[b][i], g_saxs_factors[b], coors_chunks[bb][i], g_saxs_factors[bb], binsize, nbins) >> g_radist 
    copy(g_radist) > o_radist_chunks[i]
  o_radist.join()  
  return o_radist

binsize = 0.5
nbins = 20000 #100 A * 100 A

"""
#Reducing structures for testing

xstruc = 500
xatom = 1000000
eulers = eulers[:xstruc] ###
positions = positions[:xstruc] ###
templates[0] = templates[0][:xatom] ###
templates[1] = templates[1][:xatom] ###
saxs_factors[0] = saxs_factors[0][:xatom] ###
saxs_factors[1] = saxs_factors[1][:xatom] ###
"""

def func_radist_numpy(eulers, positions, templates, saxs_factors, binsize, nbins):
  coors = func_collect(eulers, positions, templates)
  receptors, ligands = coors[0].get_data(), coors[1].get_data()

  d = receptors[:,:,None,:] - ligands[:,None,:,:]
  #e = (d*d).sum(axis=-1)
  e = np.einsum("...abc,...abc->...ab", d, d)
  w = saxs_factors[0][:,None] * saxs_factors[1][None,:]
  radist = np.array( [np.histogram(ee, bins=nbins, range=(0,nbins*binsize), weights=w)[0] for ee in e] )
  return radist


#Testing

import time
t = time.time()
ntimes = max(1, min(10, 10000/len(eulers)))
print "radist GPU, %dx, %d structures" % (ntimes, len(eulers))
for n in range(ntimes):
  radists = func_radist(eulers, positions, templates, saxs_factors, binsize, nbins).get_data()
print "Radial distribution: sum, bin index of maximum value, maximum value"
for radist in radists[:10]:
  print radist.sum(), radist.argmax(), radist[radist.argmax()]
print "radist GPU, %dx, %d structures," % (ntimes, len(eulers)), time.time() - t, "sec"

t = time.time()
maxstruc = 100
eulers = eulers[:maxstruc]
positions = positions[:maxstruc]
print "radist CPU, 1x, %d structures" % len(eulers)
radists = func_radist_numpy(eulers, positions, templates, saxs_factors, binsize, nbins)
print "Radial distribution: sum, bin index of maximum value, maximum value"
for radist in radists[:10]:
  print radist.sum(), radist.argmax(), radist[radist.argmax()]
print "radist CPU, 1x, %d structures" % len(eulers), time.time() - t, "sec"
