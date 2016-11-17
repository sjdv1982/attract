#combine_score.py

import numpy as np
import sys, os
import numpy.ma as ma
import random 
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


if '--complexes' in sys.argv:
  Benchmark = sys.argv[sys.argv.index('--complexes')+1]
  dirlist = np.genfromtxt(Benchmark,dtype=str)
  crossnum = float(sys.argv[sys.argv.index('--complexes')+2])
  if crossnum == 1.:
    crosslen = len(dirlist)
  else:
    crosslen = int(len(dirlist)*((crossnum-1.)/crossnum))
else:
  Benchmark = 'full'
  dirlist = [x for x in os.listdir('.') if os.path.isdir(x)]
  crosslen = len(dirlist)
  
numcpx=len(dirlist)

if '--method' in sys.argv:
  method=sys.argv[sys.argv.index('--method')+1]
else:
  print 'please give a metho to combine the scores'
  sys.exit()
  
if '--numstruc' in sys.argv:
  structures=int(sys.argv[sys.argv.index('--numstruc')+1])
else:
  print 'please give the number of structures to train on'
  sys.exit()
  
if '--output' in sys.argv:
  outname=sys.argv[sys.argv.index('--output')+1]
else:
  print 'please give name of the output file'
  sys.exit()

if '--scores' in sys.argv:  
  nscores=int(sys.argv[sys.argv.index('--score')+1])
  scores=sys.argv[sys.argv.index('--score')+2:sys.argv.index('--score')+2+nscores]
  head=''
  for score in scores:
    head+=score+' & '
else:
  print 'give the scores to train on'
  sys.exit()
  
counts=np.zeros((len(dirlist),structures,len(scores)))
for k,folder in enumerate(dirlist):
  for i,score in enumerate(scores):
    if score.split('.')[-1]=='rescore':
      tmp = np.genfromtxt(folder+'/'+score,usecols=[1],unpack=True)
      if len(tmp) < structures:
	print folder
	sys.exit()
	
      counts[k,:,i]=tmp[:structures]
    elif score.split('.')[-1]=='dat':
      j=0
      gobj=open(folder+'/'+score, 'r')
      glines=gobj.readlines()
      if len(glines) < structures:
	print folder
	sys.exit()
      for gline in glines:
	if 'Energy:' in gline:
	  gline=gline.strip().split()
	  counts[k,j,i]=float(gline[2])
	  j+=1
	if j==structures: break
      gobj.close()


if '--spherescores' in sys.argv:
  counts = np.sign(counts)*counts**2
  
if method=='bayes':
  from sklearn.naive_bayes import GaussianNB
  
  np.set_printoptions(precision=10,suppress=True,linewidth=2000)
  counts=counts.reshape((len(dirlist)*structures,len(scores)))
  classes=np.zeros(len(dirlist)*structures)
  y=np.zeros(len(dirlist)*structures)
  yvalue=sys.argv[sys.argv.index('--yvalue')+1]
  rcut=sys.argv[sys.argv.index('--yvalue')+2]
  kind=yvalue.split('.')[-1]
  if rcut!='multiclass':
    rcut=float(rcut)
  if rcut=='multiclass' and kind!='capstars':
    print 'Multiclass only possible with Capristars as input!!!'
    sys.exit()
  for k, folder in enumerate(dirlist):
    if kind[:5]=='irmsd' or kind[:5]=='lrmsd':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures,-1]
      y[k*structures:(k+1)*structures]=ys
      classes[k*structures:(k+1)*structures]=np.array(ys<rcut,dtype=int)
    elif kind[:4]=='fnat':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=ys
      classes[k*structures:(k+1)*structures]=np.array(ys>rcut,dtype=int)
    elif kind[:8]=='capstars':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=ys
      if rcut!='multiclass':
	classes[k*structures:(k+1)*structures]=np.array(ys<rcut,dtype=int)
      else:
	classes[k*structures:(k+1)*structures]=ys[:]

#preprocessing data for each complex alone
  from sklearn import preprocessing
  for k in range(len(dirlist)):
    scaler = preprocessing.StandardScaler().fit(counts[k*structures:(k+1)*structures])
    print scaler.mean_, scaler.std_
    counts[k*structures:(k+1)*structures]=scaler.transform(counts[k*structures:(k+1)*structures])

#for crossvalidation
  traincounts=counts[:structures*crosslen]
  trainclass=classes[:structures*crosslen]
  testcounts=counts[structures*crosslen:]
  testclass=classes[structures*crosslen:]
  
  print '##### Start fit #####'
  wclf=GaussianNB()
  wclf.fit(traincounts,trainclass)
  classprior=wclf.class_prior_
  theta=wclf.theta_
  sigma=wclf.sigma_
  #from sklearn.externals import joblib
  #joblib.dump(wclf, 'Combine_parameter_'+os.path.splitext(outname)[0]+'.pkl') 
  import cPickle
  # save the classifier
  with open('Combine_parameter_'+os.path.splitext(outname)[0]+'.pkl', 'wb') as fid:
    cPickle.dump(wclf, fid)    


  np.savetxt('Combine_parameter_'+os.path.splitext(outname)[0]+'.par', np.append(theta,sigma,axis=0) , header=head+' '+method+' '+yvalue+' '+str(rcut)+' '+Benchmark) 
  

if method=='svm':
  np.set_printoptions(precision=10,suppress=True,linewidth=2000)
  counts=counts.reshape((len(dirlist)*structures,len(scores)))
  classes=np.zeros(len(dirlist)*structures)
  y=np.zeros(len(dirlist)*structures)
  yvalue=sys.argv[sys.argv.index('--yvalue')+1]
  rcut=sys.argv[sys.argv.index('--yvalue')+2]
  kind=yvalue.split('.')[-1]
  if rcut!='multiclass':
    rcut=float(rcut)
  if rcut=='multiclass' and kind!='capstars':
    print 'Multiclass only possible with Capristars as input!!!'
    sys.exit()
  for k, folder in enumerate(dirlist):
    if kind[:5]=='irmsd' or kind[:5]=='lrmsd':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures,-1]
      y[k*structures:(k+1)*structures]=ys
      classes[k*structures:(k+1)*structures]=np.array(ys<rcut,dtype=int)
    elif kind[:4]=='fnat':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=ys
      classes[k*structures:(k+1)*structures]=np.array(ys>rcut,dtype=int)
    elif kind[:8]=='capstars':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=ys
      if rcut!='multiclass':
	classes[k*structures:(k+1)*structures]=np.array(ys<rcut,dtype=int)
      else:
	classes[k*structures:(k+1)*structures]=ys[:]
  
  if '--fitparams' in sys.argv:
    kern=sys.argv[sys.argv.index('--fitparams')+1]
    cachesize=int(sys.argv[sys.argv.index('--fitparams')+2])
    Ci=float(sys.argv[sys.argv.index('--fitparams')+3])
    prob=bool(sys.argv[sys.argv.index('--fitparams')+4])
    tolerance=float(sys.argv[sys.argv.index('--fitparams')+5]) #usually 0.001
  else:
    kern ='linear'
    cachesize=5000
    Ci=1.
    prob=False
    tolerance=0.001
#preprocessing data for each complex alone
  from sklearn import preprocessing
  for k in range(len(dirlist)):
    scaler = preprocessing.StandardScaler().fit(counts[k*structures:(k+1)*structures])
    print scaler.mean_, scaler.std_
    counts[k*structures:(k+1)*structures]=scaler.transform(counts[k*structures:(k+1)*structures])
   
  #fig=plt.figure()
  #ax=fig.add_subplot(111, projection='3d')
  #ax.scatter(counts[:,0], counts[:,1], c=classes, s=np.pi*15*(classes+0.01), alpha=0.5)
  #plt.show()
  #sys.exit()
#for crossvalidation
  traincounts=counts[:structures*crosslen]
  trainclass=classes[:structures*crosslen]
  testcounts=counts[structures*crosslen:]
  testclass=classes[structures*crosslen:]


  print '####### start SVC #########'
      
  wclf=SVC(kernel=kern, cache_size=cachesize, C=Ci, probability=prob, tol=tolerance)
  wclf.fit(traincounts,trainclass)
  if kern=='rbf' or kern=='sigmoid' or kern=='poly':
    #print wclf.support_
    supvec=wclf.support_vectors_
    print supvec
    #print wclf.n_support_
    print wclf.dual_coef_
    import cPickle
    # save the classifier
    with open('Combine_parameter_'+os.path.splitext(outname)[0]+'.pkl', 'wb') as fid:
      cPickle.dump(wclf, fid) 
    np.savetxt('Combine_parameter_'+outname, supvec , header=head+' '+method+' '+yvalue+' '+str(rcut)+' '+Benchmark+' '+kern+' probability: '+str(prob)+' tolerance: '+str(tolerance)+' C: '+str(Ci)) 

  elif kern=='linear':
    wparams=wclf.coef_
    iparams=wclf.intercept_
    np.savetxt('Combine_parameter_'+os.path.splitext(outname)[0]+'.par', -wparams, header=head+' '+method+' '+yvalue+' '+str(rcut)+' '+Benchmark)
    print -wparams
    print -iparams

    
elif method=='mc':
  #prepreparation to bring the scores all on the same dimension
  meanscores=np.ones(nscores)
  if '--meanscale' in sys.argv:	
    for s in range(nscores):
      meanscores[s]=np.mean(counts[:,:,s])
      counts[:,:,s]=counts[:,:,s]/meanscores[s]
  
  
  np.set_printoptions(precision=3,suppress=True,linewidth=2000)
  import duplicatelib
  classes=np.zeros((len(dirlist),structures))
  yvalue=sys.argv[sys.argv.index('--yvalue')+1]
  target=sys.argv[sys.argv.index('--yvalue')+2]
  y=[]
  kind=yvalue.split('.')[-1]
  for k, folder in enumerate(dirlist):
    if kind[:5]=='irmsd' or kind[:5]=='lrmsd':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures,-1]
      y.append(ys)
    elif kind[:4]=='fnat':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y.append(ys)
    elif kind[:8]=='capstars':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y.append(ys)
      
    if target=='irmsd' or target=='lrmsd':
      classes[k]=1/ys
    elif target=='fnat' or target=='capstars':
      classes[k]=ys
    else:
      tvalues=np.genfromtxt(target)
      tn=int(sys.argv[sys.argv.index('--yvalue')+3])
      for t in range(len(tvalues)-1):
	rmin=tvalues[t,0]
	rmax=tvalues[t+1,0]+0.00000001
	coor=ma.nonzero(ma.masked_outside(ys,rmin,rmax))
	for co in coor:
	  classes[k,co]=tvalues[t,tn]
	  
  tfunction=sys.argv[sys.argv.index('--fitparams')+1]
  annealing=sys.argv[sys.argv.index('--fitparams')+2]
  T0=float(sys.argv[sys.argv.index('--fitparams')+3])
  stepping=sys.argv[sys.argv.index('--fitparams')+4]
  parstep=float(sys.argv[sys.argv.index('--fitparams')+5])
  ansteps=sys.argv[sys.argv.index('--fitparams')+6]
  
  if '--duplication' in sys.argv:
    dup=int(sys.argv[sys.argv.index('--duplication')+1])
    for d,count in enumerate(counts):
      counts[d],classes[d]=duplicatelib.duplicate1bin(y[d], counts[d], target, dup)
    
  pot=20.*np.random.rand(len(scores))-10.
  print pot  
    #choice of targetfunction
  if tfunction=='refine':
      bestof=structures/10
      def targfunc(en,weight):
	evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(en, axis=-1)	#sorts position of energy after engery for every complex
	for i in range(numcpx):
	    reweight[i]=weight[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo=np.sum(reweight[:,:bestof], axis=1)/norm
	ausgabe=evaluo[:]
	return evaluo,ausgabe
  if tfunction=='rank':
      bestof=structures/10
      position=np.arange(bestof,0,-1)
      def targfunc(en,weight):
	evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(en, axis=-1)	#sorts position of energy after engery for every complex
	for i in range(numcpx):
	    reweight[i]=weight[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo=np.sum(reweight[:,:bestof]*position,axis=1)/norm
	ausgabe=np.sum(reweight[:,:bestof], axis=1)/norm
	return evaluo,ausgabe
  if tfunction=='simple':			#usually use another valuate function
      position=np.arange(1,structures+1)
      def targfunc(en,weight):
	#evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(en, axis=-1)
	for i in range(numcpx):
	  reweight[i]=weight[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo = -(ma.notmasked_edges(ma.masked_equal(reweight,0),axis=1)[0][1]+1.)
	ausgabe=evaluo[:]
	return evaluo,ausgabe
  if tfunction=='positionlinear':
      position=np.arange(structures-1,-1,-1)
      def targfunc(en,weight):
	evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(en, axis=-1)	#sorts position of energy after engery for every complex
	for i in range(numcpx):
	    reweight[i]=weight[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo=np.sum(reweight[:,:]*position,axis=1)/norm
	ausgabe=evaluo[:]
	return evaluo,ausgabe
  if tfunction=='positionquadratic':
      position=0.01*(np.arange(structures-1,-1,-1))**2
      def targfunc(en,weight):
	evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(en, axis=-1)	#sorts position of energy after engery for every complex
	for i in range(numcpx):
	    reweight[i]=weight[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo=np.sum(reweight[:,:]*position,axis=1)/norm
	ausgabe=evaluo[:]
	return evaluo,ausgabe
  if tfunction=='refine-positionlinear':
      bestof=structures/10
      position=np.arange(structures-1,-1,-1)
      def targfunc(eni,weighti):
	evaluo=np.zeros(numcpx,dtype=np.float32)
	aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
	for i in range(numcpx):
	    reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
	evaluo=(np.sum(reweight[:,:]*position,axis=1)/norm)*(np.sum(reweight[:,:bestof], axis=1)/norm)
	ausgabe=evaluo[:]
	return evaluo,ausgabe	  
	  
  #choice of annealing scheme
  if annealing=='linear':
      def function(k):
	  return (T0)-(T0/stepsfloat)*k		#linear annealing
  elif annealing=='exponential':
      def function(k):
	  return T0*(0.001)**(k/(stepsfloat))		#exponential annealing
  elif annealing=='logarithmic':
      def function(k):
	  return Emaxx/np.log(1.01+k)			#logarithmic annealing
  elif annealing=='ziczac':				#zic-zac annealing
      def function(k): 
	  return T0*(0.001)**(k/(stepsfloat))*np.cos(2*np.pi*k/(stepsfloat*0.3))**2+(T0/10.)*(0.001)**(k/(stepsfloat))
  elif annealing=='wallclimber':
      def function(k): 
	  return 1.
  #+adaptive cooling scheme coming soon
  if annealing!='wallclimber':
    def acfunc():
      return random.random()
  elif annealing=='wallclimber':
    def acfunc():
      return 1.	
	  
	  
  #choice of stepsize
  if stepping=='normal':
      def stepfunc(st,temp):
	  return random.choice((st,-st))
  elif stepping=='adaptive':
      def stepfunc(st,temp):
	  return st*random.choice((temp/T0,-temp/T0))

  en=np.zeros((numcpx,structures),dtype=np.float32)
  enold=np.zeros((numcpx,structures),dtype=np.float32)
  for j in range(numcpx):    
      for i in range(structures):
	  energy=np.sum(counts[j,i]*pot)
	  en[j,i]=energy
	  enold[j,i]=energy

  reweight=np.zeros((numcpx,structures))
  norm=np.sum(classes, axis=1)		#sums all weight up to norm the evaluation for each complex
  masknorm=ma.array(dirlist,mask=ma.getmask(ma.masked_not_equal(norm,0)))
  masknorm=ma.compressed(masknorm)
  if len(masknorm)!=0:
    print 'ATTENTION!!!: '+str(len(masknorm))+' complexes do not have any good structure to train on'
    print masknorm
    sys.exit()
    
  evaluold,ausgab=targfunc(en,classes)
  ausgeben=ausgab[:]

  crossfile=open('Combine_Annealing_'+os.path.splitext(outname)[0]+'.txt','w')
  for i in dirlist:
    crossfile.write("%7s" % i)
  crossfile.write('\n'+str(ausgeben).strip('[]'))

  #Monte-Carlo steps
  steps=int(ansteps)			#steps of Monte Carlo Annealing
  stepsfloat=float(steps)		
  Emaxx=float(len(dirlist))*0.1
  ac=0

  for k in range(steps):
      T=function(k)
      a=random.randint(0,nscores-1)
      times=stepfunc(parstep,T)	
      en[:,:]=enold+counts[:,:,a]*times
      evalu,ausgab=targfunc(en,classes)#score is the sum of the top 100 structures
      test=np.exp(-(np.sum(evaluold[:crosslen]-evalu[:crosslen]))/T)
      acceptance=acfunc()
      
      if acceptance<test:
	  pot[a]+=times
	  enold[:,:]=en
	  evaluold[:]=evalu
	  ausgeben[:]=ausgab
	  ac+=1
      #print(ausgeben)
      crossfile.write('\n'+str(ausgeben).strip('[]'))
  pot=pot/meanscores
  print pot
  np.savetxt('Combine_parameter_'+os.path.splitext(outname)[0]+'.par',pot,fmt="%6s",header=head+' '+method+' '+yvalue+' '+tfunction+' '+annealing+' '+str(T0)+' '+stepping+' '+ansteps+' '+Benchmark)	

	
elif method == 'regression':
  np.set_printoptions(precision=10,suppress=True,linewidth=2000)
  counts=counts.reshape((len(dirlist)*structures,len(scores)))
  y=np.zeros(len(dirlist)*structures)
  yvalue=sys.argv[sys.argv.index('--yvalue')+1]
  kind=yvalue.split('.')[-1]
  for k, folder in enumerate(dirlist):
    if kind[:5]=='irmsd' or kind[:5]=='lrmsd':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures,-1]
      y[k*structures:(k+1)*structures]=-1./ys
    elif kind[:4]=='fnat':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=-ys
    elif kind[:8]=='capstars':
      ys=np.genfromtxt(folder+'/'+yvalue)[:structures]
      y[k*structures:(k+1)*structures]=-ys

#preprocessing data for each complex alone
  if '--preprocessing' in sys.argv:
    prepro=sys.argv[sys.argv.index('--preprocessing')+1]
    from sklearn import preprocessing
    if prepro=='Standard':
      for k in range(len(dirlist)):
	scaler = preprocessing.StandardScaler().fit(counts[k*structures:(k+1)*structures])
	counts[k*structures:(k+1)*structures]=scaler.transform(counts[k*structures:(k+1)*structures])
    if prepro=='MinMax':
      for k in range(len(dirlist)):
	scaler = preprocessing.MinMaxScaler().fit(counts[k*structures:(k+1)*structures])
	counts[k*structures:(k+1)*structures]=scaler.transform(counts[k*structures:(k+1)*structures])
    if prepro=='complexmean':
      for k in range(len(dirlist)):
	scaler = np.mean(counts[k*structures:(k+1)*structures],axis=0)
	counts[k*structures:(k+1)*structures]=counts[k*structures:(k+1)*structures]/scaler

    
#for crossvalidation
  traincounts=counts[:structures*crosslen]
  trainclass=y[:structures*crosslen]
  testcounts=counts[structures*crosslen:]
  testclass=y[structures*crosslen:]


  print '####### start regression #########'
  
  print trainclass
  print traincounts
  
  if '--regressiontype' in sys.argv:
    regtype = sys.argv[sys.argv.index('--regressiontype')+1]
  else:
    print 'give regressiontype'
    sys.exit()

  if regtype == 'Bayesianridge':
    from sklearn.linear_model import BayesianRidge
    wclf=BayesianRidge()
    wclf.fit(traincounts,trainclass)
  
  elif regtype == 'Ridge':
    from sklearn.linear_model import Ridge
    if '--fitparams' in sys.argv:
      alp = float(sys.argv[sys.argv.index('--fitparams')+1])
    else:
      alp = 0.5
    
    wclf=Ridge(alpha=alp)			#|Wx-b|^2+alpha|W|^2
    wclf.fit(traincounts,trainclass)
  
  elif regtype == 'ols':
    from sklearn.linear_model import LinearRegression
    wclf=LinearRegression()		#regression with least squares fit
    wclf.fit(traincounts,trainclass)
  
  elif regtype=='svr-robust':
    from sklearn.svm import SVR
    if '--fitparams' in sys.argv:
      eps = float(sys.argv[sys.argv.index('--fitparams')+1])
      ci = float(sys.argv[sys.argv.index('--fitparams')+2])
    else:
      eps = 0.1
      ci = 1.
    
    wclf=SVR(kernel='linear',epsilon=eps, C=ci)  #epsilon determines the range where points can be around the fit
    wclf.fit(traincounts,trainclass)
    
    
  wparams=wclf.coef_
  iparams=wclf.intercept_  
  print wparams*10.
  np.savetxt('Combine_parameter_'+os.path.splitext(outname)[0]+'.par', wparams, header=head+' '+method+' '+yvalue+' '+regtype+' '+Benchmark)

 
  