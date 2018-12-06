import sys, os
from subprocess import Popen, PIPE
import argparse

from lib.read_struc import read_struc, depivotize
from lib.read_situs import read_situs
from lib.pdb import get_coordinates, get_coordinates_multi_attract,get_weights, get_weights_multi_attract
from lib.euler2rotmat import euler2rotmat
from lib.collect import collect
from lib.gridify import gridify
from lib.core import copy, transpose, fill, nptask
from lib.atomdensitymask import atomdensitymask1, atomdensitymask2
from axsym import AxSymAction, prepare_symtrans, apply_axsym
import numpy as np
rand = np.random.random
from topview import TVArray, TVContext, TVArrayList, TVHash

coor_chunksize = 6000
grid_chunksize = 6000 #for 20x20x20 grid
grid_chunk_parallel = 1

import sys

a = argparse.ArgumentParser(prog="neo-attract-em.py")
a.add_argument("datfile")
a.add_argument("parfile")
a.add_argument("pdbfile")
a.add_argument("pdbfile2", nargs="?")
a.add_argument("--ghost", action="store_true")
a.add_argument("--mc", action="store_true")
a.add_argument("--score", action="store_true")
a.add_argument("--mcscalerot")
a.add_argument("--mcscalecenter")
a.add_argument("--mcmax", type=int)
a.add_argument("--vmax", type=int)
a.add_argument("--mctemp", type=float, default=1.0)
a.add_argument("--gravity",type=int)
a.add_argument("--rstk", type=float)
a.add_argument("--chunks")
a.add_argument("--np")
a.add_argument("--output")
a.add_argument("--atomdensitymask", nargs=2)
a.add_argument("--rest")
a.add_argument("--axsym", nargs="*", default=[], action=AxSymAction)
args = a.parse_args()

if (not args.mc and not args.score) or not args.ghost:
  cmd = "python2 $ATTRACTDIR/../protocols/attract.py " + " ".join(sys.argv[1:])
  p = Popen(cmd, shell=True, close_fds=True)
  p.wait()
  sys.exit(0)

datfile = args.datfile
emfile = None
if args.atomdensitymask:
  emfile = args.atomdensitymask[0]
  #HACK: emfile will end with .mask, but we want the .sit file
  assert emfile.endswith(".mask")
  emfile = emfile[:-len(".mask")] + ".sit"
  assert os.path.exists(emfile)
  forceconstant = float(args.atomdensitymask[1])
pdbfiles = [v for v in (args.pdbfile, args.pdbfile2) if v is not None]
nriter = 0
if args.mc:
  nriter = int(args.mcmax)
  mcscalerot = float(args.mcscalerot)
  mcscalecenter = float(args.mcscalecenter)
  mctemp = float(args.mctemp)
rstk = args.rstk
outpf = sys.stdout
if args.output is not None:
  outpf = open(args.output, "w")


if len(pdbfiles) == 1:
  templates0 = get_coordinates_multi_attract(pdbfiles[0], ignore_weightless=False)
else:
  templates0 = [get_coordinates(pdb, ignore_weightless=False) for pdb in pdbfiles]
coms = [v.sum(axis=0)/len(v) for v in templates0]
pivots, eulers, positions0 = read_struc(datfile, coms)
positions = depivotize(pivots, eulers, positions0)
pivots = [np.array((0,0,0), dtype="float32") for p in pivots]

grid_chunksize = 6000 #for 20x20x20 grid or smaller
emdata, gridspacing, origin = None, None, None
if emfile is not None:
  emdata, gridspacing, origin = read_situs(emfile)
  origin = tuple(origin)
  s = emdata.shape[0]
  if s > 20:
    coor_chunksize = 6000
    grid_chunksize = 600

templates = templates0 #TODO: we can ignore weights in calculating the grids, but not in evaluating them...
for t,p in zip(templates, pivots):
  t -= p
if len(pdbfiles) == 1:
  weights = get_weights_multi_attract(pdbfiles[0], ignore_weightless=False)
else:
  weights = [get_weights(pdb, ignore_weightless=False) for pdb in pdbfiles]

rest, restsel = None, None
if args.rest:
  rest = read_rest(args.rest) #TODO

def overlap(eulers, positions, templates, weights, origin, gridspacing, reps_emdata, hash_emdata, coor_chunksize, grid_chunksize, maxdensity):
  """
  each element in reps_emdata: emdata multiplied by -maxdensity, replicated grid_chunksize times, and uploaded to the GPU
  """
  assert grid_chunksize <= coor_chunksize
  nbodies = len(templates)
  assert eulers.shape[1] == positions.shape[1] == nbodies
  nstruc = eulers.shape[0]
  assert eulers.shape[0] == positions.shape[0]
  assert len(reps_emdata) == grid_chunk_parallel
  for rep_emdata in reps_emdata:
    assert len(rep_emdata.shape) == 4
    assert rep_emdata.shape[0] == grid_chunksize
    assert rep_emdata.shape[1] == rep_emdata.shape[2] == rep_emdata.shape[3]

  #g_rep_emdata = TVArray("g_rep_emdata", rep_emdata, gpu=True, hash=hash_emdata)
  g_reps_emdata = TVArrayList("g_reps_emdata", reps_emdata, gpu=True, hashes=[hash_emdata]*grid_chunk_parallel)
  overlaps = TVArray("overlaps", dtype="float32", shape=(nstruc,))

  e = TVArray("e", eulers)
  p = TVArray("p", positions)
  i_templates = TVArrayList("i_templates", templates)
  i_weights = TVArrayList("i_weights", weights)
  assert len(templates) == len(weights)
  g_grids = []
  for k in range(grid_chunk_parallel):
    g_grids0 = TVArray("g_grids{%d}" % k, shape=reps_emdata[k].shape, dtype="float32", gpu=True)
    g_grids.append(g_grids0)

  g_templates = []
  g_weights = []
  for t,w in zip(i_templates, i_weights):
    tt = TVArray("g_"+t.name()[0], gpu=True)
    copy(t) > tt
    g_templates.append(tt)
    ww = TVArray("g_"+w.name()[0], gpu=True)
    copy(w) > ww
    g_weights.append(ww)

  e_chunk_T = TVArray("e_chunk_T")
  p_chunk_T = TVArray("p_chunk_T")
  g_e_chunk = TVArray("g_e_chunk", gpu=True)
  g_p_chunk = TVArray("g_p_chunk", gpu=True)
  g_rotmats = TVArray("g_rotmats", gpu=True)

  g_coors = []
  for n in range(nbodies):
    a = TVArray("g_coors{%d}" % n, gpu=True)
    g_coors.append(a)

  e_chunks = e.rchunks(coor_chunksize)
  p_chunks = p.rchunks(coor_chunksize)
  overlaps_chunks = overlaps.wchunks(coor_chunksize)
  overlaps_chunk = TVArray("overlaps_chunk", dtype="float32")

  overlaps_chunk2 = []
  for k in range(grid_chunk_parallel):
    overlaps_chunk2_0 = TVArray("overlaps_chunk2{%d}" % k, dtype="float32")
    overlaps_chunk2.append(overlaps_chunk2_0)
    #overlaps_chunk2_0.cache()

  ov, g_ov = [], []
  for k in range(grid_chunk_parallel):
    ov0 = TVArray("ov{%d}"%k)
    g_ov0 = TVArray("g_ov{%d}"%k, gpu=True)
    ov.append(ov0)
    g_ov.append(g_ov0)

  for i in range(len(e_chunks)):
    print >> sys.stderr, "CHUNK", i+1
    e_chunk, p_chunk = e_chunks[i], p_chunks[i]
    transpose(e_chunk, (1,0,2)) > e_chunk_T
    transpose(p_chunk, (1,0,2)) > p_chunk_T
    copy(e_chunk_T) > g_e_chunk
    copy(p_chunk_T) > g_p_chunk

    g_ee = g_e_chunk.rsplit()
    g_pp = g_p_chunk.rsplit()

    for n in range(nbodies):
      euler2rotmat(g_ee[n]) > g_rotmats
      collect(g_pp[n], g_rotmats, g_templates[n]) > g_coors[n]
    g_coors_chunks = [a.rchunks(grid_chunksize) for a in g_coors]

    g_e_chunk = g_e_chunk.join()
    g_p_chunk = g_p_chunk.join()

    #"overlaps_chunk" represents overlaps::i. since we cannot shard a shard

    overlaps_chunk._current().shape = (len(e_chunk),)
    overlaps_chunks2 = overlaps_chunk.wchunks(grid_chunksize)

    for j in range(len(g_coors_chunks[n])):
      k = j % grid_chunk_parallel
      overlaps_chunk2[k]._current().shape = overlaps_chunks2[j].shape
      copy(g_reps_emdata[k]) > g_grids[k]
      for n in range(nbodies):
        atomdensitymask1(g_coors_chunks[n][j], g_weights[n], origin, gridspacing) >> g_grids[k]

      fill(0) > overlaps_chunk2[k]
      for n in range(nbodies):
        atomdensitymask2(g_coors_chunks[n][j], g_grids[k], origin, gridspacing, maxdensity) > g_ov[k]
        copy(g_ov[k]) > ov[k]
        nptask(ov[k], ".sum(axis=1)") >> overlaps_chunk2[k]

      copy(overlaps_chunk2[k]) > overlaps_chunks2[j] #cannot accumulate a shard

    g_coors = [a.join() for a in g_coors]
    overlaps_chunk = overlaps_chunk.join()
    copy(overlaps_chunk) > overlaps_chunks[i] #cannot shard a shard

  overlaps = overlaps.join()
  #overlaps.cache()
  return overlaps


#prepare
import pycuda
import pycuda.autoinit
density, margin, maxdensity = None, None, None
gpu_reps_emdata, hash_emdata = None, None
if emfile:
  density = 0.824 #average protein packing density is 0.824 Da/A**3
  margin = 1.5*(gridspacing**(-2.0/3))
  maxdensity = (density+margin*density) * (gridspacing**3)
  emdata2 = np.array(emdata*-maxdensity, dtype="float32")
  hash_emdata = TVHash(emdata2)
  gpu_emdata = pycuda.gpuarray.to_gpu(emdata2)
  gpu_reps_emdata = []
  for n in range(grid_chunk_parallel):
    gpu_rep_emdata = pycuda.gpuarray.empty(dtype="float32", shape=(grid_chunksize,)+emdata.shape)
    for n in range(grid_chunksize):
      pycuda.driver.memcpy_dtod(gpu_rep_emdata[n].gpudata,gpu_emdata.gpudata, gpu_emdata.nbytes)
    gpu_reps_emdata.append(gpu_rep_emdata)
symtrans = prepare_symtrans(args.axsym, eulers.shape[1])
for snr, s in enumerate(symtrans):
  assert s.targetligand == len(templates)
  templates.append(templates[s.ligand])
  weights.append(weights[s.ligand])
pivots_axsym = [np.array((0,0,0), dtype="float32") for p in range(len(pivots) + len(symtrans))]
#mainloop

from  topview import TVContext
c_overlap  = TVContext(overlap)

def energy(eulers, positions):
  energies = np.zeros(len(eulers))
  if emfile:
    overlaps = c_overlap(eulers, positions, templates, weights, origin, \
      gridspacing, gpu_reps_emdata, hash_emdata, coor_chunksize, \
      grid_chunksize, maxdensity)
    overlaps *= forceconstant
    energies += overlaps
  if args.gravity == 1:
    gravities = (positions**2).sum(axis=2).sum(axis=1) * rstk
    energies += gravities
  return energies

def monte(eulers, positions):
  #TODO: euler2rotmat and matmult with random axis-angle matrix
  r = rand(eulers.shape) - 0.5
  r *= mcscalerot
  eulers2 = eulers + np.array(r, "float32")
  r = rand(positions.shape) - 0.5
  r *= mcscalecenter
  positions2 = positions + np.array(r, "float32")
  return eulers2, positions2

def metro(oldr, newr):
  en0, eu0, pos0 = oldr
  en1, eu1, pos1 = newr
  d = (en1 - en0).clip(min=-10, max=1000)
  r = rand(len(en1))
  d2 = np.exp(-d/mctemp) - r
  cond = (d2>0)
  en = np.where(cond, en1, en0)
  cond = cond[:, np.newaxis, np.newaxis]
  eu = np.where(cond, eu1, eu0)
  pos = np.where(cond, pos1, pos0)

  return en, eu, pos

def write_struc(f, pivots, eulers, positions, energies):
  for n in range(len(pivots)):
    print >>f,"#pivot %d %.3f %.3f %.3f" % (n+1, pivots[n][0], pivots[n][1], pivots[n][2])
  print >>f,"#centered receptor: false"
  print >>f,"#centered ligands: false"
  strucnr = 0
  for eu, pos, e in zip(eulers, positions, energies):
    strucnr += 1
    print >>f,"#" + str(strucnr)
    print >>f,"## Energy:", e
    for eu2, pos2 in zip(eu, pos):
      print >>f," ",  eu2[0], eu2[1], eu2[2], pos2[0], pos2[1], pos2[2]

eulers_axsym, positions_axsym = apply_axsym(symtrans, eulers, positions)
energies = energy(np.ascontiguousarray(eulers_axsym, dtype="float32"), np.ascontiguousarray(positions_axsym, dtype="float32"))
print >> sys.stderr, 0, sum(energies)/len(energies)
for n in range(nriter):
  eulers2, positions2 = monte(eulers, positions)
  eulers2_axsym, positions2_axsym = apply_axsym(symtrans, eulers2, positions2)
  new_energies = energy(np.ascontiguousarray(eulers2_axsym, dtype="float32"), np.ascontiguousarray(positions2_axsym, dtype="float32"))
  energies, eulers, positions = metro( (energies, eulers, positions), (new_energies, eulers2, positions2) )
  print >> sys.stderr, n+1, sum(energies)/len(energies)

write_struc(outpf, pivots, eulers, positions, energies)
