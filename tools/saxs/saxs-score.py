
from lib.read_struc import read_struc, depivotize
from lib.pdb import get_coordinates, get_vacuum_factors, get_dummy_factors
from lib.euler2rotmat import euler2rotmat
from lib.collect import collect
from lib.core import copy, transpose, fill
from lib.calc_radist import calc_radist
import numpy as np
from scipy.spatial.distance import cdist
import math
import argparse
from topview import TVArray, TVContext, TVArrayList

import sys
parser = argparse.ArgumentParser()
parser.add_argument('saxsfile',help="SAXS Data file (in Angstrom) with error column",action="store")
parser.add_argument('datfile',help="ATTRACT dat structures file",action="store")
parser.add_argument('--fit', default=False,help='Fit parameter C1',action="store_true")
parser.add_argument('--pdbs',nargs='+',required=True,help="Please provide the PDB files of your proteins in all-atom format")
args = parser.parse_args()
saxsfile = args.saxsfile
datfile = args.datfile
pdbs = args.pdbs
threads = 8
pivots, eulers, positions0 = read_struc(datfile, None)
positions = depivotize(pivots, eulers, positions0)
templates = [get_coordinates(pdb, ignore_weightless=True) for pdb in pdbs]
vacuum_factors = [get_vacuum_factors(pdb) for pdb in pdbs]
dummy_factors = [get_dummy_factors(pdb) for pdb in pdbs]

CHUNKSIZE = 2000


                   
def func_radist(eulers, positions, templates, vacuum_factors, dummy_factors, binsize, nbins):
  
  nbodies = len(templates)
  assert eulers.shape[1] == positions.shape[1] == nbodies
  nstruc = eulers.shape[0]
  assert eulers.shape[0] == positions.shape[0]
  
  o_radist = TVArray("radist", shape = (nstruc, 3, nbins), dtype="float32")
  
  e = TVArray("e", eulers)  
  p = TVArray("p", positions)
  eT = TVArray("eT")  
  pT = TVArray("pT")
  g_e = TVArray("g_e", gpu=True)  
  g_p = TVArray("g_p", gpu=True)
  g_rotmats = []
  g_coors = []
  g_templates = []
  i_templates = TVArrayList("i_templates", templates)
  g_vacuum_factors = []
  i_vacuum_factors = TVArrayList("i_vacuum_factors", vacuum_factors)  
  g_dummy_factors = []
  i_dummy_factors = TVArrayList("i_dummy_factors", dummy_factors)
  for b in range(nbodies):
    g_rotmats.append(TVArray("g_rotmats", gpu=True))
    g_coors.append(TVArray("g_coors[%d]" % b, gpu=True))
    g_templates.append(TVArray("g_templates[%d]" % b, gpu=True))
    copy(i_templates[b]) > g_templates[b]
    g_vacuum_factors.append(TVArray("g_vacuum_factors[%d]" % b, gpu=True))
    copy(i_vacuum_factors[b]) > g_vacuum_factors[b]   
    g_dummy_factors.append(TVArray("g_dummy_factors[%d]" % b, gpu=True))
    copy(i_dummy_factors[b]) > g_dummy_factors[b]
  transpose(e, (1,0,2)) > eT
  transpose(p, (1,0,2)) > pT
  copy(eT) > g_e
  copy(pT) > g_p
  g_ee = g_e.rsplit()
  g_pp = g_p.rsplit()

  for b in range(nbodies): 
    euler2rotmat(g_ee[b]) > g_rotmats[b]
    #collect(g_pp[b], g_rotmats, g_templates[b]) > g_coors[b]
      
  #coors_chunks = [g.rchunks(CHUNKSIZE) for g in g_coors]
  rotmats_chunks = [g.rchunks(CHUNKSIZE) for g in g_rotmats]
  pp_chunks = []
  for b in range(nbodies):
    g = TVArray("g_pp_body[%d]"%b, gpu=True)
    copy(g_pp[b]) > g
    pp_chunks.append(g.rchunks(CHUNKSIZE))
  o_radist_chunks = o_radist.wchunks(CHUNKSIZE)  
  g_radist = TVArray("g_radist", gpu = True, shape = (CHUNKSIZE, 3,nbins), dtype="float32")  
  for i in range(len(o_radist_chunks)):
    g_radist.shape = (rotmats_chunks[b][i].shape[0],3, nbins)
    for b in range(nbodies):
      collect(pp_chunks[b][i], rotmats_chunks[b][i], g_templates[b]) > g_coors[b]
    
    fill(0) > g_radist
    for b in range(nbodies):
      for bb in range(b+1, nbodies):
        calc_radist(g_coors[b], g_vacuum_factors[b], g_dummy_factors[b], g_coors[bb], g_vacuum_factors[bb], g_dummy_factors[bb], binsize, nbins) >> g_radist
    copy(g_radist) > o_radist_chunks[i]
  o_radist.join()  
  return o_radist
  #g_radist = TVArray("g_radist", gpu = True, shape = (nstruc, nbins), dtype="float32")
  #fill(0) > g_radist
  #for b in range(nbodies):
    #for bb in range(b+1, nbodies):
      #calc_radist(g_coors[b], g_saxs_factors[b], g_coors[bb], g_saxs_factors[bb], binsize, nbins) >> g_radist 
  #copy(g_radist) > o_radist
  #return o_radist

binsize = 0.5
nbins = 40000 #100 A * 100 A
if len(templates[0])+len(templates[1]) > 5000:
  CHUNKSIZE = 500
  threads = 4
  nbins = 80000

nbodies = len(templates)
initradist = [[],[],[]]
initradist[0] = [0.0 for jj in range(nbins)]
initradist[1] = [0.0 for jj in range(nbins)]
initradist[2] = [0.0 for jj in range(nbins)]
for n in range(nbodies):
  coors = templates[n]
  if not len(coors) == len(vacuum_factors[n]):
    print len(coors), len(vacuum_factors[n])
    sys.exit(1)
  
  Y = cdist(coors,coors,'sqeuclidean')
  w = vacuum_factors[n][:,None]*vacuum_factors[n][None,:]
  initradist[0] += np.histogram(Y, bins=nbins, range=(0,nbins*binsize),weights=w)[0]
  w = dummy_factors[n][:,None]*dummy_factors[n][None,:]
  initradist[1] += np.histogram(Y, bins=nbins, range=(0,nbins*binsize),weights=w)[0]
  w = vacuum_factors[n][:,None]*dummy_factors[n][None,:] + dummy_factors[n][:,None]*vacuum_factors[n][None,:]
  initradist[2] += np.histogram(Y, bins=nbins, range=(0,nbins*binsize),weights=w)[0]
 
bins = np.arange(0,binsize*nbins,binsize)
assert len(bins) == len(initradist[0])
q, Irefe,sigma = np.loadtxt(saxsfile,unpack=True,usecols=[0,1,2]) #TODO: read experimental profile
sigma2 = sigma**2
sincfunc = [np.sinc(qq*np.sqrt(bins)/math.pi) for qq in q]
b = 0.23 #Taken from IMP, squared b, modulation_function_parameter
E2 = np.exp(-2.0*b*(q**2)) #squared modulation function f(q) = f(0)*exp(-b*q^2)


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

def calc_chi_numpy(indata):
  radist, initradist, q, Irefe, sigma2, sincfunc, E2, C1, C1square, fit = indata
  assert len(radist) == len(initradist)
  r = 2.0*radist[0]+initradist[0]
  rr = 2.0*radist[1]+initradist[1]
  rrr = 2.0*radist[2]+initradist[2]
  chitotal = []
  chidefault = 0.0
  I1 =  [np.sum(r*sincfunc[qq]) for qq in range(len(q))]
  I2 =  [np.sum(rr*sincfunc[qq]) for qq in range(len(q))]
  I3 =  [np.sum(rrr*sincfunc[qq]) for qq in range(len(q))]
  I1 = I1*E2
  I2 =I2*E2
  I3 = I3*E2
  if fit:
    for c1, c1square in zip(C1,C1square):
      I = I1+c1square*I2-c1*I3
      #scale factor
      c = np.sum(Irefe*I/sigma2)/np.sum(I*I/sigma2)
      chi = np.sum((Irefe-I*c)**2/sigma2)
      if c1[1] - 1.0 < 0.00001:
        chidefault = chi
      chitotal = np.append(chitotal,chi)
    chifit = np.amin(chitotal)
  
  else:
    I = I1+1.0*I2-1.0*I3
    c = np.sum(Irefe*I/sigma2)/np.sum(I*I/sigma2)
    chifit = np.sum((Irefe-I*c)**2/sigma2)
    chidefault = chifit
    
  factor = 1.0/float(len(q))
  return (np.sqrt(chifit*factor), np.sqrt(chidefault*factor))
  
radists = func_radist(eulers, positions, templates, vacuum_factors, dummy_factors, binsize, nbins).get_data()
#sys.exit()
import pycuda.autoinit
close = pycuda.autoinit.context
close.detach()
nopt = 0.005
CC1 = np.arange(0.95,1.05,nopt)
C1 = []
for cc in CC1:
  tmp = [cc**3*np.exp(-(4.0/3.0*math.pi)**1.5*qq*qq*1.62*1.62*(cc**2-1.0)/(4.0*math.pi)) for qq in q]
  C1.append(np.array(tmp))
  
C1square = [cc*cc for cc in C1]
from multiprocessing import Pool
p = Pool(threads)
indata = [(radist,initradist,q,Irefe,sigma2,sincfunc,E2,C1,C1square,args.fit) for radist in radists]
results = p.map(calc_chi_numpy,indata)
p.close()
p.join()
for i,r in enumerate(results):
  print i+1, r[0], r[1]
  
