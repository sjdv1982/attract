import sys
try:
  import psyco
  psyco.full()
except:
  pass  

fast = False
try:
  i = sys.argv.index("--fast")
  sys.argv.pop(i)
  fast = True
except ValueError:
  pass  

bodies = int(sys.argv[1])
structures = int(sys.argv[2])
seed = 1
if len(sys.argv) > 3: seed = int(sys.argv[3])
import random
if seed != -1:
  random.seed(seed)

rsq = 75 * 75
from math import *

def adjust(x,y,z):
  w = sqrt(rsq/(x*x+y*y+z*z))
  return x*w,y*w,z*w

print "#pivot auto"
print "#centered receptor: true"
print "#centered ligands: true"
for n in range(structures):
  print "#"+str(n+1)
  p = []
  for nn in range(bodies):
    x = random.random()-1 
    y = random.random()-1
    z = random.random()-1
    x,y,z = adjust(x,y,z)
    p.append((x,y,z))
  if bodies == 2:
    x,y,z = p[0]
    p[1] = (-x,-y,-z)
  else:
    while 1:  
      delta = 0
      newp = []
      for i in range(bodies):
	x,y,z = p[i]
	xold,yold,zold = x,y,z
	for j in range(bodies):
          if i==j: continue
          xx,yy,zz = p[j]
	  dx,dy,dz = xold-xx,yold-yy,zold-zz
	  dsq = dx*dx+dy*dy+dz*dz
	  dsq2 = 30.0/dsq
	  x += dx * dsq2
	  y += dy * dsq2
	  z += dz * dsq2
	  dx,dy,dz = x-xx,y-yy,z-zz
	  #dsqnew = dx*dx+dy*dy+dz*dz
	  #print dsq,dsqnew-dsq
	x,y,z = adjust(x,y,z)
	if not fast:
          dd = x-xold,y-yold,z-zold
	  delta+=dd[0]*dd[0]+dd[1]*dd[1]+dd[2]*dd[2]
	newp.append((x,y,z))
      p = newp
      #print delta
      if fast or delta < 0.001*bodies: break
  for i in range(bodies):
    x,y,z = p[i]
    print 2*(random.random()-1)*pi,2*(random.random()-1)*pi,2*(random.random()-1)*pi,x,y,z
