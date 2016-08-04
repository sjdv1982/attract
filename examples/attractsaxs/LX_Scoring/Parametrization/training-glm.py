import numpy as np
from sklearn.linear_model import LinearRegression
import sys, os
import numpy.ma as ma
import duplicatelib

Header = ''
if '--complexes' in sys.argv:
  Complexes=sys.argv[sys.argv.index('--complexes')+1]	#set of complexes to train on
  crossnum = int(sys.argv[sys.argv.index('--complexes')+2])
  Header += ' Benchmark: '+Complexes+' crossnum: '+str(crossnum)
else:
  print 'Please give a list of complexes to train on'
  sys.exit()
  
if '--grid' in sys.argv:
  Counts = sys.argv[sys.argv.index('--grid')+1]	#precalculated grids for the structures of the complexes
  strux = sys.argv[sys.argv.index('--grid')+2]
  gridtype = sys.argv[sys.argv.index('--grid')+3]
  Header += ' Grid: '+Counts+' #structures: '+strux+' gridtype: '+gridtype
  mctype = 'normal'
  if gridtype == 'distances':
    mctype = sys.argv[sys.argv.index('--grid')+4]
    rbincut = int(sys.argv[sys.argv.index('--grid')+5])
    bincut = int(sys.argv[sys.argv.index('--grid')+6])
    binlist = np.zeros(bincut, dtype=int) 
    Header += ' rcutbin: '+str(rbincut)+' bins: '
    for i in range(bincut):
      binlist[i] = int(sys.argv[sys.argv.index('--grid')+7+i])-1
      Header += ' '+str(binlist[i])
else:
  print 'please give a grid to train on'
  sys.exit()
  
  
if '--regressiontype' in sys.argv:
  regressiontype = sys.argv[sys.argv.index('--regressiontype')+1]
else:
  print 'insert type of regression. give parameter for fitting with --fitparams' 
  sys.exit()
#Robustlinearmodels
#svr-robust
#ols
#nonneglsq
#Ridge
#BayesianRidge
#Lasso
#elasticNet
#RANSAC
#logistic
#SGDClass

prepro = False
if '--preprocessing' in sys.argv:
  prepro = True
  preprotype = sys.argv[sys.argv.index('--preprocessing')+1]
else: 
  print 'preprocessing of complexes with very different number of contatcs could be important !!!!'
  
if '--evaluate' in sys.argv:
  Rmsd=sys.argv[sys.argv.index('--evaluate')+1]
  dup=sys.argv[sys.argv.index('--evaluate')+2]
  Header += ' weight: '+dup
  if dup == 'duplication' or dup == 'probabilities':
    valuate=sys.argv[sys.argv.index('--evaluate')+3]
    prorow = int(sys.argv[sys.argv.index('--evaluate')+4])
    Header += valuate+' '+str(prorow)
else:
  print 'insert file to evaluate structures on'
  sys.exit()
  
if '--cutoffeval' in sys.argv:
  evcut = float(sys.argv[sys.argv.index('--cutoffeval')+1])

if '--maxweight' in sys.argv:
  maxcut = float(sys.argv[sys.argv.index('--maxweight')+1])
    
if '--functionshape' in sys.argv:
  bined = int(sys.argv[sys.argv.index('--functionshape')+1])
  atyps = int(sys.argv[sys.argv.index('--functionshape')+2])
  if bined == 0:
    bins = 1
    partypes = atyps
  else:
    bins = bined
    partypes=(atyps*(atyps-1))/2+atyps
else:
  print 'give number of bins and number of atomtypes'
  sys.exit()
  
if '--eraseatomtype' in sys.argv:
  eratomtype = int(sys.argv[sys.argv.index('--eraseatomtype')+1])
  hydrogens = True
  Header += ' Eraseatomtype: '+str(eratomtype)
else:
  hydrogens = False

if '--insertnatives' in sys.argv:
  natives=True
  nativegrid=sys.argv[sys.argv.index('--insertnatives')+1]
  Header += ' insert natives '
else:
  natives=False

if '--output' in sys.argv:
  output = sys.argv[sys.argv.index('--output')+1]
  outform = sys.argv[sys.argv.index('--output')+2]
else: 
  print 'insert name for outputfiles'
  sys.exit()

if outform == 'saddle':
  if '--powers' in sys.argv:
    power=np.ones(bincut)
    for i in range(bincut):
      power[i] = float(sys.argv[sys.argv.index('--powers')+1+i])
  else:
    print 'for saddle output one has to give the power of the bins'
    sys.exit()

weighton = False
weighting = 0
if '--weighting' in sys.argv:
  weighton = True
  weighting=int(sys.argv[sys.argv.index('--weighting')+1])

predicton = False
if '--predict' in sys.argv:
  predicton = True

#insert posibility to put some features nonlinearly together like [12,13]-> z529=x12*x13

cmplx=np.genfromtxt(Complexes, dtype=str)
numcpx=len(cmplx)
structures=int(strux)

if crossnum == 1:
  trainlen = numcpx*structures
else:
  trainlen = int(((crossnum-1.)*numcpx)/crossnum)*structures
testlen = numcpx-trainlen

rmsd=np.zeros((numcpx,structures),dtype=np.float32)
#read in rmsd values
for d,name in enumerate(cmplx):					#for every complex one makes a matrix with all the counts of atompairs
    #read rmsds
    fobj=np.genfromtxt(name+'/'+Rmsd)[:structures]
    if len(np.shape(fobj))>1:
      rmsd[d*structures:(d+1)*structures]=fobj[:structures,1]
    else:
      rmsd[d*structures:(d+1)*structures]=fobj[:structures]   
    if natives:
      if Rmsd.split('.')[-1][1:5]=='rmsd':
	rmsd[d,-1]=0.001
      elif Rmsd.split('.')[-1][:8]=='capstars':
	rmsd[d,-1]==4
      elif Rmsd.split('.')[-1][:4]=='fnat':
	rmsd[d,-1]==1.
      else:
	print 'rmsd type for natives not understood'
	sys.exit()

allcounts = np.zeros((structures*numcpx, bins*partypes),dtype = np.float32)
weight=np.zeros((numcpx*structures), dtype = np.float32)

for d, name in enumerate(cmplx):
  if gridtype == 'distances':
    distcounts = np.load(name+'/'+Counts, mmap_mode='r')
    distcounts = distcounts[:rbincut,:,:structures,:]      
    counts=np.sum(distcounts, axis=0)[binlist]
    if mctype == 'keepsign':
      for i in range(bincut):
	if (i+1)%2 == 0:
	  counts[i]=-counts[i]
  elif gridtype == 'step':
    counts = np.load(name+'/'+Counts, mmap_mode='r')
    if len(np.shape(counts)) < 3:
      counts = counts[:structures,:]
      counts = counts.reshape((1,structures,partypes))
    elif len(np.shape(counts)) == 3:
      counts = counts[:,:structures,:]
    else:
      print 'shape of grids not understood', np.shape(counts)
      sys.exit()
  
  if natives:
    if gridtype == 'distances':
      natdistcounts = np.load(name+'/'+nativegrid, mmap_mode='r')
      natdistcounts = distcounts[:rbincut]
      natcounts=np.sum(distcounts, axis=0)[binlist]
      if mctype == 'keepsign':
	for i in range(bincut):
	  if (i+1)%2 == 0:
	    natcounts[i]=-natcounts[i]
    elif gridtype == 'step':
      natcounts = np.load(name+'/'+nativegrid, mmap_mode='r')
      if len(np.shape(natcounts)) < 3:
	natcounts = natcounts.reshape((1,structures,partypes))
    counts=np.append(counts[:,:-1,:],natcounts, axis=1)
    
  if dup == 'duplication':
    counts,weight[d*structures:(d+1)*structures]=duplicatelib.duplicate(rmsd[d,:], counts, valuate, prorow)
  for z in range(bins):
    allcounts[d*structures:(d+1)*structures, z*partypes: (z+1)*partypes] = counts[z]

rmsd = rmsd.reshape(numcpx*structures)
if dup == 'irmsd' or dup == 'lrmsd':
  weight[:]=1./rmsd[:]
  if '--cutoffweight' in sys.argv:
    evcut=1./evcut
  if '--maxweight' in sys.argv:
    maxcut = 1./maxcut
elif dup=='fnat' or dup=='capstars':
  weight[:]=rmsd[:]
elif dup == 'probabilities':
  targetvalues=np.genfromtxt(valuate)	#weights for targetfunction
  for i in range(len(targetvalues)-1):
    rmin=targetvalues[i,0]
    rmax=targetvalues[i+1,0]
    maske=ma.masked_outside(rmsd,rmin,rmax+0.00001)
    coor=ma.nonzero(maske)
    for j in range(len(coor[0])):
      weight[coor[0][j]]=targetvalues[i,prorow]		#targetvalues are in the file 1-3 prob, 4 weight

if '--cutoffweight' in sys.argv:
    weight=ma.masked_less(weight,evcut).filled(0.)

if '--maxweight' in sys.argv:
  weight = ma.masked_greater(weight, maxcut).filled(maxcut)

if prepro:
    from sklearn import preprocessing
    if preprotype=='Standard': # for each complex it scales over each parameter (axis = 0)!!!!
      for d in range(numcpx):
	scaler=preprocessing.StandardScaler().fit(allcounts[d*structures:(d+1)*structures])
	allcounts[d*structures:(d+1)*structures]=scaler.transform(allcounts[d*structures:(d+1)*structures])
    elif preprotype=='MinMax': #for every structure is makes the features between 0,1
      for d in range(numcpx):
	scaler=preprocessing.MinMaxScaler().fit(allcounts[d*structures:(d+1)*structures])
	allcounts[d*structures:(d+1)*structures]=scaler.transform(allcounts[d*structures:(d+1)*structures])
    elif preprotype=='meancomplex':
      for d in range(numcpx): #divides contacts by mean sum of all contacts for a structures in a complex 
	scaler = np.mean(np.sum(allcounts[d*structures:(d+1)*structures], axis = 1))
	allcounts[d*structures:(d+1)*structures]=allcounts[d*structures:(d+1)*structures]/scaler
    else:
      print '#### no scaler chosen for preprocessing ####'
      sys.exit()

if '--eraseatomtype' in sys.argv:
    hydindex=0
    hydindexes=[]
    for n in range(atyps):
	for nn in range(n,atyps):
	    if nn==eratomtype-1:        
	      for ap in range(bins):
		hydindexes.append(ap*partypes+hydindex)
	    hydindex+=1
    allcounts=np.delete(counts, hydindexes, -1)
    partypes=np.shape(traincounts)[-1]
    atyps =atyps - 1

#Regression:
weight = -weight
trainweight = weight[:trainlen]
traincounts = allcounts[:trainlen]
testweight = weight[trainlen:]
testcounts = allcounts[trainlen:]



############### all generalized linear models for fitting ################

if regressiontype=='Robustlinearmodels':
    import statsmodels.api as sm
    
    traincounts2=np.append(traincounts,np.ones((len(traincounts),1)),axis=1)

    if '--regressionparams' in sys.argv:
      Mweight = sys.argv[sys.argv.index('--regressionparams')+1]
      if Mweight == 'Huber':
	Mweight = sm.robust.norms.HuberT()
      elif Mweight == 'Andrew':
	Mweight = sm.robust.norms.AndrewWave()
      elif Mweight == 'Hampel':
	Mweight = sm.robust.norms.Hampel()
      elif Mweight == 'Ramsay':
	Mweight = sm.robust.norms.RamsayE()
      elif Mweight == 'TrimmedMean':
	Mweight = sm.robust.norms.TrimmedMean()
      elif Mweight == 'Tukey':
	Mweight = sm.robust.norms.TukeyBiweight()
    else:
      Mweight = sm.robust.norms.HuberT()
    
    wclf = sm.RLM(trainweight, traincounts2, M=Mweight)
    
    if '--fitparams' in sys.argv:
      conv = sys.argv[sys.argv.index('--fitparams')+1] #coefs, weights, sresid, dev
      if conv == 'help':
	help(sm.RLM.fit)
	sys.exit()
      cov = sys.argv[sys.argv.index('--fitparams')+2] #H1, H2, H3, H4
      maxiter = int(sys.argv[sys.argv.index('--fitparams')+3])
      scale_est = sys.argv[sys.argv.index('--fitparams')+4] #mad, stand_mad, Huber
      if scale_est == 'Huber':
	scales_est = sm.robust.scale.HuberScale()
      tol = float(sys.argv[sys.argv.index('--fitparams')+5])
      rlm_results = wclf.fit(conv = conv, cov = cov, maxiter = maxiter, scale_est = scale_est, tol = tol)
    
    else:
      rlm_results = wclf.fit()
    
    allparams=rlm_results.params
    wparams=allparams[:bins*partypes]
    iparams=allparams[bins*partypes]



elif regressiontype=='svr-robust':
    from sklearn.svm import SVR
    
    if '--regressionparams' in sys.argv:
      eps = float(sys.argv[sys.argv.index('--regressiontype')+1])
      ci = float(sys.argv[sys.argv.index('--regressiontype')+2])
      tolerance = float(sys.argv[sys.argv.index('--regressiontype')+3])
      cache_size = float(sys.argv[sys.argv.index('--regressiontype')+4])
      max_iter = int(sys.argv[sys.argv.index('--regressiontype')+5])
    else: 
      eps = 0.01
      ci = 10000.
      tolerance = 0.001
      cache_size = 4000
      max_iter = -1
      
    wclf=SVR(kernel='linear',epsilon=eps, C=ci, tol = tolerance, cache_size = cache_size, max_inter = max_iter)  #epsilon determines the range where points can be around the fit
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_


elif regressiontype=='ols':
    
    traincounts2=np.append(traincounts,np.ones((len(traincounts),1)),axis=1)
    
    allparams=np.linalg.lstsq(traincounts2,trainweight)[0]
    wparams=allparams[:bins*partypes]
    iparams=allparams[bins*partypes]

    
elif regressiontype=='nonneglsq':
    
    print 'start.........'

    from scipy.optimize import nnls

    traincounts2=np.append(traincounts,np.ones((len(traincounts),1)),axis=1)
    allparams,r=nnls(traincounts2,trainweight) 
    wparams=allparams[:bins*partypes]
    iparams=allparams[bins*partypes]
    

elif regressiontype=='Ridge':
    
    if '--regressionparams' in sys.argv:
      alph = float(sys.argv[sys.argv.index('--regressiontype')+1])
      solve = sys.argv[sys.argv.index('--regressiontype')+2] #auto, svd, cholesky, lsqr, sparse_cg
      tolerance = float(sys.argv[sys.argv.index('--regressiontype')+3])
      max_iter = int(sys.argv[sys.argv.index('--regressiontype')+4])
    else: 
      alph = 1.
      solve = 'auto'
      tolerance = 0.001
      max_iter = None
      
    from sklearn import linear_model
    wclf=linear_model.Ridge(alpha=.5)			#|Wx-b|^2+alpha|W|^2
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_



elif regressiontype=='Lasso':
  
    if '--regressionparams' in sys.argv:
      alph = float(sys.argv[sys.argv.index('--regressiontype')+1])
      positive = bool(sys.argv[sys.argv.index('--regressiontype')+2])
      tolerance = float(sys.argv[sys.argv.index('--regressiontype')+3])
      max_iter = int(sys.argv[sys.argv.index('--regressiontype')+4])
    else: 
      alph = 1.
      positive = False
      tolerance = 0.0001
      max_iter = 1000
  
    from sklearn import linear_model
    wclf=linear_model.Lasso(alpha=alph, positive = positive, max_iter = max_iter, tol = tolerance)			#(1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_    

elif regressiontype=='elasticNet':
      
    if '--regressionparams' in sys.argv:
      alph = float(sys.argv[sys.argv.index('--regressiontype')+1])
      positive = bool(sys.argv[sys.argv.index('--regressiontype')+2])
      tolerance = float(sys.argv[sys.argv.index('--regressiontype')+3])
      max_iter = int(sys.argv[sys.argv.index('--regressiontype')+4])
      l1_ratio = int(sys.argv[sys.argv.index('--regressiontype')+5])
    else: 
      alph = 1.
      positive = False
      tolerance = 0.0001
      max_iter = 1000
      l1_ratio = 0.5
    
    from sklearn import linear_model
    wclf=linear_model.ElasticNet(alpha = alph, tol = tolerance, max_iter = max_iter, l1_ratio = l1_ratio, positive = positive)		#1 / (2 * n_samples) * ||y - Xw||^2_2 + alpha * l1_ratio * ||w||_1+ 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_

elif regressiontype=='RANSAC':
    from sklearn.linear_model import RANSACRegressor			#good to use if outliers are in the dataset
    
    if '--regressionparams' in sys.argv:
      #base_estimator=None
     # min_samples=None
      #residual_threshold=None
      #is_data_valid=None
      #is_model_valid=None
      max_trials=int(sys.argv[sys.argv.index('--regressionparams')+1])
      #stop_n_inliers=inf
      #stop_score=inf
      stop_probability=float(sys.argv[sys.argv.index('--regressionparams')+2])
      #residual_metric=None
      #random_state=None
      wclf = RANSACRegressor( max_trials = max_trials, stop_probability = stop_probability)
    else:
      wclf = RANSACRegressor()
    
    wclf.fit(traincounts,trainweight)
    wparams=wclf.estimator_.coef_
    iparams=wclf.estimator_.intercept_

elif regressiontype=='BayesianRidge':
    from sklearn.linear_model import BayesianRidge
    
    if '--regressionparams' in sys.argv:
      n_iter=int(sys.argv[sys.argv.index('--regressiontype')+1])
      tol=float(sys.argv[sys.argv.index('--regressiontype')+2])
      alpha_1=float(sys.argv[sys.argv.index('--regressiontype')+3])
      alpha_2=float(sys.argv[sys.argv.index('--regressiontype')+4])
      lambda_1=float(sys.argv[sys.argv.index('--regressiontype')+5])
      lambda_2=float(sys.argv[sys.argv.index('--regressiontype')+6])
      normalize=bool(sys.argv[sys.argv.index('--regressiontype')+7])
      wclf=BayesianRidge(n_inter = n_inter, tol = tol, alpha_1 = alpha_1, alpha_2 = alpha_2, lambda_1 = lambda_1, lambda_2 = lambda_2, normalize = normalize)
    else: 
      wclf=BayesianRidge()
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_


if regressiontype=='logistic':
    from sklearn import linear_model
    classborder=float(sys.argv[sys.argv.index('--regressiontype')+2])

    if dup == 'irmsd' or dup == 'lrmsd':
      classborder=-1./classborder
    elif dup=='fnat' or dup=='capstars':
      classborder=-classborder    
    trainweight=np.array(trainweight<=classborder)
    
    if '--regressionparams' in sys.argv:
      classweight=int(sys.argv[sys.argv.index('--fitparams')+1])
      Ci=float(sys.argv[sys.argv.index('--fitparams')+2])
      tolerance=float(sys.argv[sys.argv.index('--fitparams')+3]) #usually 0.001
      penalty = sys.argv[sys.argv.index('--fitparams')+4]  #l1 or l2
      wclf=linear_model.LogisticRegression(penalty = penalty, C = Ci, tol = tolerance, class_weight={1:classweight} )
    else:
      wclf=linear_model.LogisticRegression() #same algorithm as SVC(linear), LinearSVC()
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_


elif regressiontype=='SGDClass':
    from sklearn.linear_model import SGDClassifier
    classborder=float(sys.argv[sys.argv.index('--regressiontype')+2])

    if dup == 'irmsd' or dup == 'lrmsd':
      classborder=-1./classborder
    elif dup=='fnat' or dup=='capstars':
      classborder=-classborder    
    trainweight=np.array(trainweight<=classborder)
    if '--regressionparams' in sys.argv:
      loss = sys.argv[sys.argv.index('--fitparams')+1]
      #loss : str, hinge, log, modified_huber, squared_hinge, perceptron 
      # or a regression loss: squared_loss, huber, epsilon_insensitive, or squared_epsilon_insensitive 
      #The loss function to be used. Defaults to hinge, which gives a linear SVM.
      #The log loss gives logistic regression, a probabilistic classifier.
      #modified_huber is another smooth loss that brings tolerance to outliers as well as probability estimates. 
      #squared_hinge is like hinge but is quadratically penalized. perceptron is the linear loss used by the perceptron algorithm.
      #The other losses are designed for regression but can be useful in classification as well; see SGDRegressor for a description
      alpha = float(sys.argv[sys.argv.index('--fitparams')+2])
      l1_ratio = float(sys.argv[sys.argv.index('--fitparams')+3])
      classweight=int(sys.argv[sys.argv.index('--fitparams')+4])  
      penalty = sys.argv[sys.argv.index('--fitparams')+5]  
      epsilon = float(sys.argv[sys.argv.index('--fitparams')+6])
      wclf = SGDClassifier(loss = loss, alpha = alpha, l1_ratio = l1_raito, class_weight={1:classweight}, penalty = penalty, epsilon = epsilon )
    else:
      wclf = SGDClassifier()  #Gradient descent method for SVM/ probably the same as SVC()
    #wclf.transform(median) #can also reduce the number of features to the most important
    wclf.fit(traincounts,trainweight)
    wparams=wclf.coef_
    iparams=wclf.intercept_

    
########### GIVE THE PARAMETER OUT ################  

if bined == 0:
  paramatrix = wparams
  print paramatrix
  ending = ''
else:
  if outform == 'saddle':
    paramatrix = np.ones(((bins+1)*atyps,atyps))
    outpar = np.zeros((bins, partypes))
    ABpar = wparams.reshape((bins, partypes))
    outpar[0] = (ABpar[0]/ABpar[1])**(1./(power[0]-power[1]))
    outpar[1] = (ABpar[0]**(-power[1]/(power[0]-power[1])))*(ABpar[1]**(power[0]/(power[0]-power[1])))
    ending = '_parm'
  else:
    ending = ''
    if mctype == 'keepsign':
      outpar = wparams.reshape((bins,partypes))
      paramatrix=np.zeros((bins*atyps,atyps))

      for nnn in bins:
	if (nnn+1)%2 == 0:
	  outpar[nnn]=-outpar[nnn]

    else:
      outpar = wparams.reshape((bins,partypes))
      paramatrix=np.zeros((bins*atyps,atyps))
    
    ind = 0
    for n in range(atyps):
      for nn in range(n,atyps):
	for nnn in range(bins):
	  paramatrix[nnn*atyps+n,nn]=outpar[nnn,ind]
	  paramatrix[nnn*atyps+nn,n]=outpar[nnn,ind]
	ind+=1      
    
paramatrix=np.around(paramatrix,decimals=8)
np.savetxt('LinReg-Parameter_'+output+ending+'.par', paramatrix, fmt="%10s", header='Linear Regression:'+Header)    
    
        
if predicton:
      prediction=np.sum(wparams*testcounts,axis=-1)+iparams
      predictiontrain=np.sum(wparams*traincounts,axis=-1)+iparams
      outtrain=np.append(predictiontrain.reshape((len(predictiontrain),1)),trainweight.reshape((len(predictiontrain),1)),axis=1)
      out=np.append(prediction.reshape((len(prediction),1)),testweight.reshape((len(prediction),1)),axis=1)
      np.savetxt('LinReg-Predictions_testset_'+output+'.dat',out,fmt="%6s")
      np.savetxt('LinReg-Predictions_training_'+output+'.dat',outtrain,fmt="%6s")
    
