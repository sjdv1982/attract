import numpy as np
import sys
q = []
I = []
error = []
if len(sys.argv) > 2:
  if not (sys.argv[2] == '--noerror' or sys.argv[2] == '--convert'):
    print "Option not recognized!"
    sys.exit(1)
    
  if sys.argv[2] == '--noerror':
    q,I  = np.loadtxt(sys.argv[1],unpack=True)
    error = np.zeros(len(q))
    
  else:
    q,I,error = np.loadtxt(sys.argv[1],unpack=True)
  
else:
  q,I,error = np.loadtxt(sys.argv[1],unpack=True)

if len(sys.argv)> 2 and sys.argv[2] == '--convert':
  q = q/10.0
  
randnr = np.abs(np.random.poisson(10,len(q))/10.0-1.0)+1.0
weight = randnr*I*0.03*5.0*(q+0.001)
#print np.mean(error/weight)
newq, newI, newerror = [], [], []
for i in range(len(q)):
  if 5*weight[i] < error[i]:
    ##if error[i] > 0.05*I[i]:
    continue
    ##print q[i], I[i],error[i]
   
  if weight[i] < error[i]:
    weight[i] = error[i]
  newq.append(q[i])
  newI.append(I[i])
  newerror.append(weight[i])
  
#datapoints = int(sys.argv[2])
#chunks = len(newq)/datapoints
#if chunks == 0:
  #for i,j,k in zip(newq,newI,newerror):
    #print i,j,k
    
  #sys.exit()
  
#newq = newq[::chunks]
#newI = newI[::chunks]
#newerror = newerror[::chunks]
for i,j,k in zip(newq,newI,newerror):
  print i,j,k