import numpy as np
import numpy.ma as ma

def duplicate(rmsd,grid,prob,tar):
      #targetfunction
    targetvalues=np.genfromtxt(prob)	#weights for targetfunction
    weight=np.zeros(len(rmsd))
    for i in range(len(targetvalues)-1):
	rmin=targetvalues[i,0]
	rmax=targetvalues[i+1,0]
	maske=ma.masked_outside(rmsd,rmin,rmax)
	coor=ma.nonzero(maske)
	for j in range(len(coor)):
	    weight[coor[j]]=targetvalues[i,tar]
    pos=np.arange(len(rmsd))
    norm=np.sum(weight)
    duplications=ma.masked_equal(weight,0)
    dupmask=ma.getmask(duplications)
    position=ma.compressed(ma.array(pos,mask=dupmask))
    proba=ma.compressed(duplications)/norm
    duplates=(1./10.)*len(rmsd)
    n=proba*duplates
    an=np.ceil(n)
    weighting=n/an
#    print proba, n, an, weighting
    ind=0
    counts=np.zeros(np.shape(grid))
    newweight=np.zeros(len(rmsd))
    for z,i in enumerate(an):
      for j in range(int(i)):
	counts[:,ind,:]=grid[:,position[z],:]
	newweight[ind]=weighting[z]
	ind+=1
    fillpos=ma.compressed(ma.array(pos,mask=~dupmask))
    ind2=0
#    print ind
#    print np.shape(fillpos)
    for z in range(ind,len(rmsd)):
      counts[:,z,:]=grid[:,fillpos[ind2],:]
      ind2+=1
    return counts, newweight
  
def duplicate1bin(rmsd,grid,prob,tar):
      #targetfunction
    targetvalues=np.genfromtxt(prob)	#weights for targetfunction
    weight=np.zeros(len(rmsd))
    for i in range(len(targetvalues)-1):
	rmin=targetvalues[i,0]
	rmax=targetvalues[i+1,0]
	maske=ma.masked_outside(rmsd,rmin,rmax)
	coor=ma.nonzero(maske)
	for j in range(len(coor)):
	    weight[coor[j]]=targetvalues[i,tar]
    pos=np.arange(len(rmsd))
    norm=np.sum(weight)
    duplications=ma.masked_equal(weight,0)
    dupmask=ma.getmask(duplications)
    position=ma.compressed(ma.array(pos,mask=dupmask))
    proba=ma.compressed(duplications)/norm
    duplates=(1./10.)*len(rmsd)
    n=proba*duplates
    an=np.ceil(n)
    weighting=n/an
#    print proba, n, an, weighting
    ind=0
    counts=np.zeros(np.shape(grid))
    newweight=np.zeros(len(rmsd))
    for z,i in enumerate(an):
      for j in range(int(i)):
	counts[ind,:]=grid[position[z],:]
	newweight[ind]=weighting[z]
	ind+=1
    fillpos=ma.compressed(ma.array(pos,mask=~dupmask))
    ind2=0
#    print ind
#    print np.shape(fillpos)
    for z in range(ind,len(rmsd)):
      counts[z,:]=grid[fillpos[ind2],:]
      ind2+=1
    return counts, newweight
