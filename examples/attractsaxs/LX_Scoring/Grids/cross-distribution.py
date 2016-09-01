import sys, os
import numpy as np

if '--complexes' in sys.argv:
  complexes = np.genfromtxt(sys.argv[sys.argv.index('--complexes')+1], dtype = str)
else:
  complexes = [x for x in os.listdir('.') if os.path.isdir(x)]
  complexes = np.array(complexes)
  
if '--numtestset' in sys.argv:
  numtestset=int(sys.argv[sys.argv.index('--numtestset')+1])
else:
  print 'insert number of testsets for crossvalidation'
  sys.exit()
  
if '--exceptions' in sys.argv:
  numexept = int(sys.argv[sys.argv.index('--exceptions')+1])
  excepts=sys.argv[sys.argv.index('--exceptions')+2:sys.argv.index('--exceptions')+2+numexept]
  for excep in excepts:
    for i, check in enumerate(complexes):
      if check==excep:
	complexes = np.delete(complexes, i)

if '--output' in sys.argv:
  outname=sys.argv[sys.argv.index('--output')+1]
else:
  print 'please give outputfile a name'
  sys.exit()
  
cmplen=len(complexes)

permut = np.arange(cmplen)
permut = [permut[j-int(cmplen/numtestset)] for j in range(cmplen)]

for i in range(numtestset):
    complexes=complexes[permut]
    print complexes
    np.savetxt(outname+'-cv'+str(i+1)+'.txt',complexes, fmt="%s")
