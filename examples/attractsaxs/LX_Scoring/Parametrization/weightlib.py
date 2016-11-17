#weightlib.py
from scipy.spatial.distance import cdist, pdist
import numpy as np


def weightfunc(counts, weights, weightingfunc):
    if len(counts)<=25000:
	if weightingfunc==-1:
	    xdist=cdist(counts, counts, 'euclidean')+np.eye(len(counts))
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights,weights, 'euclidean')
	    ones=np.ones(np.shape(ydist))
	    onesvec=np.ones(np.shape(weights))
	    cd=ones-ydist
	    norm=np.sum(cd,axis=0)
	    ce= cd/xdist
	    cf= np.sum(ce,axis=0)
	    cg=cf/norm
	    a=abs(weights).reshape((1,len(counts)))
	    weight=cg*a
	elif weightingfunc==-2:
	    a=np.exp(-(1.-abs(weights))**2)
	    weight=a.reshape((1,len(counts)))
	elif weightingfunc==-3:
	    xdist=cdist(counts, counts, 'euclidean')+np.eye(len(counts))
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights,weights, 'euclidean')
	    ones=np.ones(np.shape(ydist))
	    onesvec=np.ones(np.shape(weights))
	    cd=ones-ydist
	    norm=np.sum(cd,axis=0)
	    ce= cd/xdist
	    cf= np.sum(ce,axis=0)
	    cg=cf/norm
	    weight=cg.reshape((1,len(cg)))
	elif weightingfunc==-4:
	    strucs=[(100*i) for i in range(len(counts)/100)]
	    xdist=xdist=cdist(counts[strucs], counts, 'euclidean')
	    sort=np.sort(xdist, axis=1)
	    sigma=np.mean(sort[:,len(counts)/1000],axis=0)
	    print sigma
	    xdist=cdist(counts, counts, 'sqeuclidean')
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights,weights, 'euclidean')
	    ones=np.ones(np.shape(ydist))
	    onesvec=np.ones(np.shape(weights))
	    cd=ones-ydist
	    norm=np.sum(cd,axis=0)
	    gauss=np.exp(-(xdist)/sigma**2)-np.eye(len(counts))
	    ce= cd*gauss
	    cf= np.sum(ce,axis=0)
	    cg=cf/norm
	    weight=cg.reshape((1,len(cg)))
	elif weightingfunc==-5:
	    strucs=[(100*i) for i in range(len(counts)/100)]
	    xdist=xdist=cdist(counts[strucs], counts, 'euclidean')
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights[strucs],weights, 'euclidean')
	    sortx=np.median(xdist, axis=1)
	    sigma=np.median(sortx,axis=0)
	    sorty=np.median(ydist, axis=1)
	    sigmay=np.median(sorty,axis=0)
	    print sigma
	    print sigmay
	    lamda=0.5
	    print lamda
	    xdist=cdist(counts, counts, 'sqeuclidean')
	    ydist=cdist(weights,weights, 'sqeuclidean')
	    cd=np.exp(-ydist/sigmay**2)
	    norm=np.sum(cd,axis=0)
	    gauss=np.exp(-(xdist)/sigma**2)-lamda*np.eye(len(counts))
	    ce= cd*gauss
	    cf= np.sum(ce,axis=0)
	    cg=cf/norm
	    weight=cg.reshape((1,len(cg)))
	elif weightingfunc==-6:
	    import math
	    def student(n,x):
		return (math.gamma(n+1./2.)/np.sqrt(n*np.pi)*math.gamma(n/2.))*(1+(x**2/n))**(-(n+1.)/2)
	    xdist=cdist(counts, counts, 'euclidean')
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights,weights, 'euclidean')
	    cd=student(70,ydist)
	    norm=np.sum(cd,axis=1)
	    gauss=student(50,xdist)
	    ce= cd*gauss
	    cf= np.sum(ce,axis=1)
	    cg=cf/norm
	    weight=cg.reshape((1,len(cg)))
	elif weightingfunc>0:
	    weight=weightingfunc*np.ones((1,len(counts)))
    else:
	if weightingfunc==-1:
	    weight=np.zeros((1,len(counts)))
	    for i in range(len(counts)/100):
		xdist=cdist(counts[i*100:(i+1)*100], counts, 'euclidean')+np.eye(len(counts))[i*100:(i+1)*100]
		weights=weights.reshape(len(weights),1)
		ydist=cdist(weights[i*100:(i+1)*100],weights, 'euclidean')
		ones=np.ones(np.shape(ydist))
		cd=ones-ydist
		norm=np.sum(cd,axis=1)
		ce= cd/xdist
		cf= np.sum(ce,axis=1)
		cg=cf/norm
		weight[0,i*100:(i+1)*100]=cg.reshape((1,len(cg)))
	    a=abs(weights).reshape((1,len(counts)))
	    weight=cg*a
	elif weightingfunc==-2:
	    a=np.exp(-(1.-abs(weights))**2)
	    weight=a.reshape((1,len(counts)))
	elif weightingfunc==-3:
	    weight=np.zeros((1,len(counts)))
	    for i in range(len(counts)/100):
		xdist=cdist(counts[i*100:(i+1)*100], counts, 'euclidean')+np.eye(len(counts))[i*100:(i+1)*100]
		weights=weights.reshape(len(weights),1)
		ydist=cdist(weights[i*100:(i+1)*100],weights, 'euclidean')
		ones=np.ones(np.shape(ydist))
		cd=ones-ydist
		norm=np.sum(cd,axis=1)
		ce= cd/xdist
		cf= np.sum(ce,axis=1)
		cg=cf/norm
		weight[0,i*100:(i+1)*100]=cg.reshape((1,len(cg)))
	elif weightingfunc==-4:
	    weight=np.zeros((1,len(counts)))
	    strucs=[(100*i) for i in range(len(counts)/100)]
	    xdist=xdist=cdist(counts[strucs], counts, 'euclidean')
	    sort=np.sort(xdist, axis=1)
	    sigma=np.mean(sort[:,len(counts)/1000],axis=0)
	    print sigma
	    for i in range(len(counts)/100):
		xdist=cdist(counts[i*100:(i+1)*100], counts, 'sqeuclidean')
		weights=weights.reshape(len(weights),1)
		ydist=cdist(weights[i*100:(i+1)*100],weights, 'euclidean')
		ones=np.ones(np.shape(ydist))
		cd=ones-ydist
		norm=np.sum(cd,axis=1)
		ident=np.zeros(np.shape(xdist))
		ident[:,i*100:(i+1)*100]=np.eye(100)
		gauss=np.exp(-(xdist)/sigma**2)-ident
		ce= cd*gauss
		cf= np.sum(ce,axis=1)
		cg=cf/norm
		weight[0,i*100:(i+1)*100]=cg.reshape((1,len(cg)))
	elif weightingfunc==-5:
	    weight=np.zeros((1,len(counts)))
	    strucs=[(100*i) for i in range(len(counts)/100)]
	    xdist=xdist=cdist(counts[strucs], counts, 'euclidean')
	    weights=weights.reshape(len(weights),1)
	    ydist=cdist(weights[strucs],weights, 'euclidean')
	    sortx=np.median(xdist, axis=1)
	    sigma=np.median(sortx,axis=0)
	    sorty=np.median(ydist, axis=1)
	    sigmay=np.median(sorty,axis=0)
	    print sigma
	    print sigmay
	    lamda=0.5
	    print lamda
	    for i in range(len(counts)/100):
		xdist=cdist(counts[i*100:(i+1)*100], counts, 'sqeuclidean')
		ydist=cdist(weights[i*100:(i+1)*100],weights, 'sqeuclidean')
		cd=np.exp(-(ydist/sigmay**2))
		norm=np.sum(cd,axis=1)
		ident=np.zeros(np.shape(xdist))
		ident[:,i*100:(i+1)*100]=np.eye(100)
		gauss=np.exp(-(xdist)/sigma**2)-lamda*ident
		ce= cd*gauss
		cf= np.sum(ce,axis=1)
		cg=cf/norm
		weight[0,i*100:(i+1)*100]=cg.reshape((1,len(cg)))
	elif weightingfunc==-6:
	    import math
	    weight=np.zeros((1,len(counts)))
	    def student(n,x):
		return (math.gamma(n+1./2.)/np.sqrt(n*np.pi)*math.gamma(n/2.))*(1+(x**2/n))**(-(n+1.)/2)
	    for i in range(len(counts)/10):
		xdist=cdist(counts[i*10:(i+1)*10], counts, 'euclidean')
		weights=weights.reshape(len(weights),1)
		ydist=cdist(weights[i*10:(i+1)*10],weights, 'euclidean')
		cd=student(70,ydist)
		norm=np.sum(cd,axis=1)
		gauss=student(50,xdist)
		ce= cd*gauss
		cf= np.sum(ce,axis=1)
		cg=cf/norm
		weight[0,i*10:(i+1)*10]=cg.reshape((1,len(cg)))
	elif weightingfunc>0:
	    weight=weightingfunc*np.ones((1,len(counts)))
    return weight