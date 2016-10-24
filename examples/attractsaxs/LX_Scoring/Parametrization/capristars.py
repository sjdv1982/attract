#capristars.py

import numpy as np
import sys, os
import numpy.ma as ma


rmsd = np.genfromtxt(sys.argv[1])[:,-1]
irmsd = np.genfromtxt(sys.argv[2])[:,-1]
fnat = np.genfromtxt(sys.argv[3])

starvec=np.zeros(np.shape(fnat),dtype=np.int8)
acceptfnat=ma.getmask(ma.masked_outside(fnat,0.09999,0.3))
mediumfnat=ma.getmask(ma.masked_outside(fnat,0.29999,0.5))
highfnat=ma.getmask(ma.masked_less_equal(fnat,0.5))
acceptfnatrmsd=ma.array(rmsd,mask=acceptfnat)
acceptfnatirmsd=ma.array(irmsd,mask=acceptfnat)
onestarrmsd=ma.nonzero(ma.masked_greater(acceptfnatrmsd, 10.0))
onestarirmsd=ma.nonzero(ma.masked_greater(acceptfnatirmsd, 4.0))
for i in onestarrmsd:
  starvec[i]=1
for i in onestarirmsd:
  starvec[i]=1
mediumfnatrmsd=ma.array(rmsd,mask=mediumfnat)
mediumfnatirmsd=ma.array(irmsd,mask=mediumfnat)
stillone=ma.nonzero(mediumfnatrmsd)
for i in stillone:
  starvec[i]=1
twostarrmsd=ma.nonzero(ma.masked_greater(acceptfnatrmsd, 5.0))
twostarirmsd=ma.nonzero(ma.masked_greater(acceptfnatirmsd, 2.0))
for i in twostarrmsd:
  starvec[i]=2
for i in twostarirmsd:
  starvec[i]=2
highfnatrmsd=ma.array(rmsd,mask=highfnat)
highfnatirmsd=ma.array(irmsd,mask=highfnat)
stilltwo=ma.nonzero(highfnatrmsd)
for i in stilltwo:
  starvec[i]=2
threestarrmsd=ma.nonzero(ma.masked_greater(highfnatrmsd, 1.0))
threestarirmsd=ma.nonzero(ma.masked_greater(highfnatirmsd, 1.0))
for i in threestarrmsd:
  starvec[i]=3
for i in threestarirmsd:
  starvec[i]=3
np.savetxt(os.path.splitext(sys.argv[1])[0]+".capstars",starvec,fmt="%d")