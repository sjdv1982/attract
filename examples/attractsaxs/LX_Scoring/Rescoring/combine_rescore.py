#combine_rescore.py
import numpy as np
import sys,os
from sklearn import preprocessing

method=sys.argv[1]
combinepar=sys.argv[2]

if '--output' in sys.argv:
  outname=os.path.splitext(sys.argv[sys.argv.index('--output')+1])[0]
else:
  outname=os.path.splitext(combinepar)[0]
  

if '--scores' in sys.argv:
    print '!!! other scores taken !!!'
    nscores = int(sys.argv[sys.argv.index('--scores')+1])
    line = sys.argv[sys.argv.index('--scores')+2:sys.argv.index('--scores')+2+nscores]
else: 
    obj=open(combinepar,'r')
    line=obj.readline()
    line=line.strip('# ').split(' & ')
    line=line[:-1]

  
scorecounts=[]
for r in range(len(line)):
    score=line[r]

    if score.split('.')[-1]=='rescore':
      scorecounts.append(np.genfromtxt(score))
    elif score.split('.')[-1]=='dat':
      gobj=open(score, 'r')
      glines=gobj.readlines()
      scorecount=[]
      for gline in glines:
	if gline[3:10]=='Energy:':
	  gline=gline.strip().split()
	  scorecount.append(float(gline[2]))
	gobj.close()
      scorecounts.append(scorecount)
scorecounts=np.array(scorecounts)
scorecounts=scorecounts.T

if '--natives' in sys.argv:  
  natscorecounts=[]
  for r in range(len(line)):
      score=line[r]
      natscorecounts.append(np.genfromtxt(os.path.splitext(score)[0]+'-native.rescore'))
  natscorecounts=np.array(natscorecounts)
  natscorecounts=natscorecounts.T

if '--spherescores' in sys.argv:
  scorecounts = np.sign(scorecounts)*scorecounts**2

if '--preprocessing' in sys.arv:
  premode=sys.argv[sys.argv.index('--preprocessing')+1]
  if premode=='Standard':
    scaler=preprocessing.StandardScaler().fit(scorecounts)
    scorecounts=scaler.transform(scorecounts)
    if '--natives' in sys.argv:
      natscorecounts=scaler.transform(natscorecounts)
  elif premode=='MinMax':
    scaler=preprocessing.MinMaxScaler().fit(scorecounts)
    scorecounts=scaler.transform(scorecounts)
    if '--natives' in sys.argv:
      natscorecounts=scaler.transform(natscorecounts)
  elif premode=='complexmean':
    scaler=np.mean(scorecounts,axis=0)
    scorecounts=scorecounts/scaler
    if '--natives' in sys.argv:
      natscorecounts=natscorecounts/scaler
  else:
    print 'no preprocessing mode given'
    sys.exit()
    

if method=='linear':
  
  if combinepar == 'average':
    compar=np.ones(len(line))
  else:
    compar=np.genfromtxt(combinepar, skip_header=1)
  
  if '--natives' in sys.argv:
    natnewscore=np.sum(natscorecounts*compar) 
  else:
    newscore=np.sum(scorecounts*compar,axis=1)

  
elif method=='nonlinear':
  scaler=preprocessing.StandardScaler().fit(scorecounts)
  scorecounts=scaler.transform(scorecounts)
  # load it again
  import cPickle
  with open(combinepar, 'rb') as fid:
    wclf = cPickle.load(fid)
    
  if '--probability' in sys.argv:
    newscore=wclf.predict_proba(scorecounts)[:,0]
  else:
    newscore=-wclf.decision_function(scorecounts)

if '--natives' in sys.argv:
  print os.path.splitext(outname)[0]+'-native.rescore'
  np.savetxt(os.path.splitext(outname)[0]+'-native.rescore',[natnewscore],fmt="%4.6f")  
else:
  print os.path.splitext(outname)[0]+'.rescore'
  np.savetxt(os.path.splitext(outname)[0]+'.rescore',newscore,fmt="%4.6f")

