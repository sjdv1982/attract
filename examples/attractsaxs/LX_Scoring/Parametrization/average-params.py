#average-params.py

import sys, os
import numpy as np

if '--scores' in sys.argv:
  ind = sys.argv.index('--scores')
  numcross = int(sys.argv[ind+1])
else: 
  print 'please give scores to average'
  sys.exit()
  
if '--bins' in sys.argv:
  bins = int(sys.argv[sys.argv.index('--bins')+1])
else:
  print 'please give number of bins, 0 for 1d array'
  sys.exit()

if '--output' in sys.argv:
  outname = sys.argv[sys.argv.index('--output')+1]
else:
  print 'please insert name for output file'
  sys.exit()

scores=[]
scoresname=''
for i in range(numcross):
  scores.append(np.genfromtxt(sys.argv[ind +2 +i], skip_header=1))
  scoresname+=' & '+sys.argv[ind +2 +i]
scores = np.array(scores)


if '--scaleparams' in sys.argv:
  if bins == 0:
    std = np.std(scores, axis = 1)
    for i in range(numcross):
      scores[i] = scores[i]/std[i]
  else: 
    std = np.ones((numcross, bins))
    atyps = len(scores[0,0])
    for i in range(numcross):
      for j in range(bins):
	scores[i,j*atyps:(j+1)*atyps]= scores[i,j*atyps:(j+1)*atyps]/np.std(scores[i,j*atyps:(j+1)*atyps])
  outparams = np.mean(scores, axis = 0)

else:

  outparams = np.mean(scores, axis = 0)

print outparams

np.savetxt(outname, outparams, header = scoresname)

