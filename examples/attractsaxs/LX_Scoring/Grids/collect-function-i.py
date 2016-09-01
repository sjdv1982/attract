import sys, os
from scipy.spatial.distance import cdist
import numpy as np
import numpy.ma as ma
import collectlibpy as collectlib
import collectgridlib as cg

#run it by: python collect-tobicount.py out_unbound-sorted-dr.dat out_unbound.npy ubAtob.pdb ubBtob.pdb
dat=sys.argv[1]

if '--proteinmodel' in sys.argv:
  atomtype=sys.argv[sys.argv.index('--proteinmodel')+1]
  receptor=sys.argv[sys.argv.index('--proteinmodel')+2]
  ligand=sys.argv[sys.argv.index('--proteinmodel')+3]
else:
  print 'please insert a proteinmodel'
  sys.exit()
  
if '--gridtype' in sys.argv:
  polation=sys.argv[sys.argv.index('--gridtype')+1]
else:
  print 'please insert gridtype'
  sys.exit()
  
if '--output' in sys.argv:
  countfile = sys.argv[sys.argv.index('--output')+1]
else:
  countfile ='pregrid_'+atomtype+'_'+os.path.splitext(dat)[0]+'_'+polation+'.npy'

single = True
if '--double' in sys.argv:
  single = False


collectinsert=[dat,receptor,ligand]

modefile=None
name=None
ensfiles=[]
anr = 0
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
    
    if anr <= len(sys.argv)-2 and arg == "--name":
      name = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue

if modefile: collectinsert+=['--modes', modefile]
#if namefile: collectinsert+=['--name', namefile]

for nr, ensfile in ensfiles:
  collectinsert += ["--ens", nr, ensfile]

collectlib.collect_init(collectinsert)

fobj=open(receptor, 'r')
lines=fobj.readlines()
linelen=len(lines)
recattyp=[]
for i,line in enumerate(lines):
  line=line.strip()
  if line[:4]=='ATOM':
    recattyp.append(int(line[57:59]))
fobj.close()

fobj=open(ligand, 'r')
lines=fobj.readlines()
linelen=len(lines)
ligattyp=[]
for i,line in enumerate(lines):
  line=line.strip()
  if line[:4]=='ATOM':
    ligattyp.append(int(line[57:59]))
fobj.close()
    
if atomtype=='tobi':
    mp = np.arange(0,33)
    replace={32:19}
    mp[replace.keys()] = replace.values()
    ligattyp = mp[ligattyp]
    recattyp = mp[recattyp]
    parcount = 19
elif atomtype=='opls':
    mp = np.arange(0,83)
    replace={30:2, 31:3, 65:4, 66:5, 67:6, 68:7, 69:8, 70:9, 71:10, 80:11, 81:12, 82:13}
    mp[replace.keys()] = replace.values()
    ligattyp = mp[ligattyp]
    recattyp = mp[recattyp]
    parcount = 13
elif atomtype=='attract':   
    mp = np.arange(0,100)
    replace={99:32}
    mp[replace.keys()] = replace.values()
    ligattyp = mp[ligattyp]
    recattyp = mp[recattyp]
    parcount = 32
elif atomtype == 'gaa':
    parcount = 27
elif atomtype == 'undefined':
    parcount = int(sys.argv[sys.argv.index('--proteinmodel')+4])
else:
  print 'proteinmodel not understood'
  sys.exit()

indexmapping = {}
ind = 0
for n in range(parcount):
  for nn in range(n,parcount):
    indexmapping[n+1,nn+1] = ind
    indexmapping[nn+1,n+1] = ind
    ind += 1

if '--maxstruc' in sys.argv:
  maxstruc=int(sys.argv[sys.argv.index('--maxstruc')+1])
else:
  maxstruc=100000
  
nparams = ind  
rlen = len(recattyp)
nstruc=0

if polation=='spline':			#interpolates so that x is always at the left border, except from the last few points
    # spline interpolation with lagrange basis, ninterpol = 1 means linear spline interpolation, 2 quadratic ...
    if '--polationparams' in sys.argv:
      interpol=sys.argv[sys.argv.index('--polationparams')+1]
      if interpol == 'none':
	ninterpol = 1
      elif interpol == 'linear':
	ninterpol = 2
      elif interpol == 'quadratic':
	ninterpol = 3 
      elif interpol == 'cubic':
	ninterpol = 4
      else:
	ninterpol = int(interpol)
      stepsize=float(sys.argv[sys.argv.index('--polationparams')+2])
      start=float(sys.argv[sys.argv.index('--polationparams')+3]) #shifts contacts automatically higher than start parameter if less then start
      end=float(sys.argv[sys.argv.index('--polationparams')+4])
    else:
      ninterpol = 4
      stepsize = 0.35
      start = 1.4
      end = 7.
    
    potsteps = int((end-start)/stepsize)+1
    ai = False
    if '--attractpar' in sys.argv:
      import fgenlib
      ai = True
      attpar = np.genfromtxt(sys.argv[sys.argv.index('--attractpar')+1], skip_header=1)
      lenpar = len(attpar[0])
      obj = open(sys.argv[sys.argv.index('--attractpar')+1], 'r')
      power1 = float(obj.readline()[:2])
      obj.close()
      eps = np.zeros((parcount, parcount))
      sig = np.zeros((parcount, parcount))
      ivor = np.ones((parcount, parcount))
      sig[:lenpar,:lenpar] = attpar[:lenpar]
      eps[:lenpar,:lenpar] = attpar[lenpar:2*lenpar]
      ivor[:lenpar,:lenpar] = attpar[lenpar*2:3*lenpar]
      parlen = parcount*(parcount-1)/2+parcount
      epsi = np.zeros(parlen)
      sigi = np.zeros(parlen)
      ivori = np.zeros(parlen)
      ind = 0
      for n in range(parcount):
	for nn in range(n,parcount):
	  epsi[ind] = eps[n,nn]
	  sigi[ind] = sig[n,nn]
	  ivori[ind] = ivor[n,nn]
	  ind+=1
      atpars = np.array([sigi, epsi, ivori])
      par=fgenlib.fgen(stepsize, start, end, potsteps, 0, atpars, power1, 6, nparams, -1, 0, 0) 

    counts = np.zeros( shape=(potsteps, maxstruc, nparams))
    while 1:
        if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
          
	result = collectlib.collect_next()
	if result or nstruc == maxstruc: break
	nstruc+=1
	coor = collectlib.collect_all_coor()
	coor_receptor = coor[:rlen]
	coor_ligand = coor[rlen:]
	counts[:,nstruc-1,:] = cg.spline(ninterpol, stepsize, start, end, potsteps, coor_receptor, coor_ligand, recattyp, ligattyp, nparams)
	if ai:
	  print np.sum(counts[:,nstruc-1,:]*par)
	  if nstruc == maxstruc: sys.exit()
    counts = counts[:,:nstruc,:] #to reduce the matrix from 100000 to nstruc

if polation=='none':
    if '--polationparams' in sys.argv:
      stepsize=float(sys.argv[sys.argv.index('--polationparams')+1])
      start=float(sys.argv[sys.argv.index('--polationparams')+2])#shifts contacts automatically higher than start parameter if less then start
      end=float(sys.argv[sys.argv.index('--polationparams')+3])
      potsteps = int((end-start)/stepsize)
      def gridcount(rcoor, lcoor):
	grid=cg.gridnone(stepsize, start, end,potsteps, rcoor, lcoor, recattyp, ligattyp, nparams)
	return grid
    elif '--bins' in sys.argv:
      numbin=int(sys.argv[sys.argv.index('--bins')+1])
      cutoffs=np.zeros(numbin)
      potsteps=numbin
      for b in range(numbin):
	cutoffs[b]=float(sys.argv[sys.argv.index('--bins')+2+b])**2
      def gridcount(rcoor, lcoor):
	grid=cg.stepgrid(cutoffs, rcoor, lcoor, recattyp, ligattyp, nparams)
	return grid
    else:
      print 'please give whether --polationparms or --bins to define the ranges of the steps'
      sys.exit()
    
    counts = np.zeros(dtype=np.uint32, shape=(potsteps, maxstruc, nparams))
    while 1:
        if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	
	result = collectlib.collect_next()
	if result or nstruc==maxstruc: break
	nstruc+=1
	coor = collectlib.collect_all_coor()
	coor_receptor = coor[:rlen]
	coor_ligand = coor[rlen:]
	counts[:,nstruc-1,:]=gridcount(coor_receptor, coor_ligand)
    counts = counts[:,:nstruc,:] #to reduce the matrix from 100000 to nstruc

if polation=='distances':
  # for opls use --powers 2 12 6
    if '--polationparams' in sys.argv:
      stepsize=float(sys.argv[sys.argv.index('--polationparams')+1])
      start=float(sys.argv[sys.argv.index('--polationparams')+2])#shifts contacts automatically higher than start parameter if less then start
      end=float(sys.argv[sys.argv.index('--polationparams')+3])
    else:
      stepsize = 1.
      start = 2.
      end = 10.
    potsteps = int((end-start)/stepsize)
    if '--powers' in sys.argv:
      numpower = int(sys.argv[sys.argv.index('--powers')+1])
      power = np.zeros(numpower)
      for i in range(numpower):
	power[i]=float(sys.argv[sys.argv.index('--powers')+2+i])
    else:
      numpower = 6
      power=np.arange(7.,1.,-1.)*2.
    counts = np.zeros(shape=(potsteps,numpower, maxstruc, nparams))
    power = power/2.
    while 1:
        if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
          
	result = collectlib.collect_next()
	if result or nstruc == maxstruc: break
	nstruc+=1
	coor = collectlib.collect_all_coor()
	coor_receptor = coor[:rlen]
	coor_ligand = coor[rlen:]
	counts[:,:,nstruc-1,:]=cg.distances(power, stepsize, start, end, potsteps, coor_receptor, coor_ligand, recattyp, ligattyp, nparams)
      
    counts = counts[:,:,:nstruc,:] #to reduce the matrix from 100000 to nstruc
	
if single and polation!='none':	
  counts = counts.astype(np.float32)

np.save(countfile, counts)

