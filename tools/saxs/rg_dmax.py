import rmsdlib
import sys
import numpy, math
from scipy.spatial.distance import cdist

anr = 0
output = None

while 1:
  anr += 1
      
  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]
  
  if anr <= len(sys.argv)-2 and arg == "--output":
    output = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
  
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')

count =1  
for unbounds in rmsdlib.read_multi_pdb(sys.argv[1]):
  coor = []  
  for p in unbounds: 
    p.remove_hydrogens()
    for c in p.coordinates():
      coor.append(c)
  coor = numpy.array(coor)
  Y2 = cdist(coor,coor,'sqeuclidean')
  rgyr2 = numpy.sum(Y2)
  rgyr2 /=len(coor)
  rgyr2 /=len(coor)
  rgyr2 *= 0.5
  dmax = math.sqrt(numpy.amax(Y2))
  print count,math.sqrt(rgyr2), dmax
  count +=1
  
  