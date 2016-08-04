import sys
import numpy as np
import fgenlib

counting=sys.argv[1]
params=sys.argv[2]

if '--gridtype' in sys.argv:
  gridtype = sys.argv[sys.argv.index('--gridtype')+1]
  bins = int(sys.argv[sys.argv.index('--gridtype')+2])
  attyps = int(sys.argv[sys.argv.index('--gridtype')+3])
  partyps = (attyps*(attyps-1))/2+attyps
  if gridtype == 'distances':
    rcut = int(sys.argv[sys.argv.index('--gridtype')+4])
    bincut = int(sys.argv[sys.argv.index('--gridtype')+5])    
    binlist = np.zeros(bincut, dtype=int)
    for i in range(bincut):
      binlist[i] = int(sys.argv[sys.argv.index('--gridtype')+6+i])-1
  elif gridtype == 'interpolate':
    import fgenlib
    start = float(sys.argv[sys.argv.index('--gridtype')+4])
    rcutoff = float(sys.argv[sys.argv.index('--gridtype')+5])
    stepsize = float(sys.argv[sys.argv.index('--gridtype')+6])
    fgentype = int(sys.argv[sys.argv.index('--gridtype')+7])
    power1 = float(sys.argv[sys.argv.index('--gridtype')+8])
    power2 = float(sys.argv[sys.argv.index('--gridtype')+9])
    potsteps = int((rcutoff-start)/stepsize)+1
else:
  print 'please give information on the gridtype'
  sys.exit()


#Tobiparams=sys.argv[3]

counts=np.load(counting,mmap_mode='r')

#tobiparams=np.genfromtxt(Tobiparams,skip_header=1)
if gridtype=='nonlinear':
  import cPickle
  with open(params.replace('txt', 'pkl'), 'rb') as fid:
    wclf = cPickle.load(fid)
elif gridtype == 'interpolate':
  parameter = np.genfromtxt(params, skip_header = 1)
  params = np.zeros((bins, partyps))
  ind=0
  for n in range(attyps):
    for nn in range(n,attyps):
      for b in range(bins):
	params[b,ind]=parameter[b*attyps+n,nn]
      ind+=1
  par = pot = fgenlib.fgen(stepsize, start, rcutoff, potsteps, fgentype, params, power1, power2, partyps, -1, 1,0.2)
else:
  parameter=np.genfromtxt(params,skip_header=1)
  if bins == 0:
    par = parameter[:]
  else:
    par=np.zeros((bins,partyps))
    ind=0
    for n in range(attyps):
      for nn in range(n,attyps):
	for b in range(bins):
	  par[b,ind]=parameter[b*attyps+n,nn]
	ind+=1

if gridtype == 'distances':
  counts = counts[:rcut]
  counts = np.sum(counts,axis = 0)[binlist]

if gridtype=='nonlinear':
  scorecounts=np.zeros((np.shape(counts)[-2],bins*partyps))
  for z in range(bins):
    scorecounts[:,z*partyps:(z+1)*partyps]=counts[z,:,:]

  if '--preprocessing' in sys.argv:
    from sklearn import preprocessing
    if sys.argv[sys.argv.index('--preprocessing')+1]=='Standard':
      scaler=preprocessing.StandardScaler().fit(scorecounts)
      scorecounts=scaler.transform(scorecounts)
    elif sys.argv[sys.argv.index('--preprocessing')+1]=='MinMax':
      scaler=preprocessing.MinMaxScaler().fit(scorecounts)
      scorecounts=scaler.transform(scorecounts)
    else:
      print '#### no scaler chosen for preprocessing ####'
      sys.exit()

  if '--probability' in sys.argv:
    newscore=wclf.predict_proba(scorecounts)[:,0]
  else:
    newscore=-wclf.decision_function(scorecounts)
  for i in range(len(scorecounts)):     
      print "%10s" % round(newscore[i],5)

else:
  if bins == 0:
    struc=len(counts)
    for i in range(struc):
      score=np.sum(counts[i,:]*par)
      score=round(score,5)
    # tobiscore=np.sum(counts[:,i,:]*tobipar)	
      print "%10s" % score
  else:
    struc=len(counts[0])
    for i in range(struc):
      score=np.sum(counts[:,i,:]*par)
      score=round(score,5)
    # tobiscore=np.sum(counts[:,i,:]*tobipar)	
      print "%10s" % score

