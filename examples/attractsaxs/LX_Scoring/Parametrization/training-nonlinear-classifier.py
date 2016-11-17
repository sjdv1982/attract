import numpy as np
from sklearn.svm import SVC
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
  if gridtype == 'distances':
    rbincut = int(sys.argv[sys.argv.index('--grid')+4])
    bincut = int(sys.argv[sys.argv.index('--grid')+5])
    binlist = np.zeros(bincut, dtype=int) 
    Header += ' rcutbin: '+str(rbincut)+' bins: '
    for i in range(bincut):
      binlist[i] = int(sys.argv[sys.argv.index('--grid')+6+i])-1
      Header += ' '+str(binlist[i])
else:
  print 'please give a grid to train on'
  sys.exit()

  
if '--classifier' in sys.argv:
  classificationtype = sys.argv[sys.argv.index('--classifier')+1]
else:
  print 'insert type of regression. give parameter for fitting with --fitparams' 
  sys.exit()

prepro = False
if '--preprocessing' in sys.argv:
  prepro = True
  preprotype = sys.argv[sys.argv.index('--preprocessing')+1]
else:
  print 'grids have to be preprocessed before one can train on them'
  sys.exit()

  
if '--evaluate' in sys.argv:
  Rmsd=sys.argv[sys.argv.index('--evaluate')+1]
  rmsdcut = float(sys.argv[sys.argv.index('--evaluate')+2])
  dup=sys.argv[sys.argv.index('--evaluate')+3]
  Header += ' weight: '+dup
  if dup == 'duplication' or dup == 'probabilities':
    valuate=sys.argv[sys.argv.index('--evaluate')+4]
    prorow = int(sys.argv[sys.argv.index('--evaluate')+5])
    Header += valuate+' '+str(prorow)
else:
  print 'insert file to evaluate structures on'
  sys.exit()
  
if '--cutoffeval' in sys.argv:
  evcut = float(sys.argv[sys.argv.index('--cutoffeval')+1])

if '--maxweight' in sys.argv:
  maxcut = float(sys.argv[sys.argv.index('--maxweight')+1])
    
if '--functionshape' in sys.argv:
  bins = int(sys.argv[sys.argv.index('--functionshape')+1])
  atyps = int(sys.argv[sys.argv.index('--functionshape')+2])
  if bins == 0:
    partypes = atyps
  else:
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

weighton = False
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
  elif gridtype == 'step':
    counts = np.load(name+'/'+Counts, mmap_mode='r')
    if len(np.shape(counts)) < 3:
      counts = counts[:structures,:]
      counts = counts.reshape((bins,structures,paratyps))
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
    elif gridtype == 'step':
      natcounts = np.load(name+'/'+nativegrid, mmap_mode='r')
      if len(np.shape(natcounts)) < 3:
	natcounts = natcounts.reshape((bins,structures,paratyps))
    counts=np.append(counts[:,:-1,:],natcounts, axis=1)
    
  if dup == 'duplication':
    counts,weight[d*structures:(d+1)*structures]=duplicatelib.duplicate(rmsd[d,:], counts, valuate, prorow)
  for z in range(bins):
    allcounts[d*structures:(d+1)*structures, z*partypes: (z+1)*partypes] = counts[z]

rmsd = rmsd.reshape(numcpx*structures)
if dup == 'irmsd' or dup == 'lrmsd':
  weight[:]=1./rmsd[:]
  rmsdcut = 1./rmsdcut
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
      scaler=preprocessing.MinMaxScaler().fit(allcounts)
      allcounts=scaler.transform(allcounts)
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
    paratyps=np.shape(traincounts)[-1]
    atyps =atyps - 1

#Regression:
weight = np.array(weight > rmsdcut, dtype = int)
trainweight = weight[:trainlen]
traincounts = allcounts[:trainlen]
testweight = weight[trainlen:]
testcounts = allcounts[trainlen:]


if classificationtype=='logistic':
    from sklearn import linear_model
    
    if '--fitparams' in sys.argv:
      classweight=int(sys.argv[sys.argv.index('--fitparams')+1])
      Ci=float(sys.argv[sys.argv.index('--fitparams')+2])
      tolerance=float(sys.argv[sys.argv.index('--fitparams')+3]) #usually 0.001
      penalty = sys.argv[sys.argv.index('--fitparams')+4]  #l1 or l2
      wclf=linear_model.LogisticRegression(penalty = penalty, C = Ci, tol = tolerance, class_weight={1:classweight} )
    else:
      wclf=linear_model.LogisticRegression() #same algorithm as SVC(linear), LinearSVC()



elif classificationtype=='SGDClass':
    from sklearn.linear_model import SGDClassifier
    
    if '--fitparams' in sys.argv:
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



elif classificationtype == 'SVM':
  
  if '--fitparams' in sys.argv:
    classweight=int(sys.argv[sys.argv.index('--fitparams')+1])
    kernelfunc=sys.argv[sys.argv.index('--fitparams')+2]
    probabilities=bool(sys.argv[sys.argv.index('--fitparams')+3])
    miter=int(sys.argv[sys.argv.index('--fitparams')+4])
    cachesize=int(sys.argv[sys.argv.index('--fitparams')+5])
    Ci=float(sys.argv[sys.argv.index('--fitparams')+6])
    tolerance=float(sys.argv[sys.argv.index('--fitparams')+7]) #usually 0.001
    
    wclf=SVC(kernel=kernelfunc,class_weight={1:classweight},probability=probabilities,max_iter=miter, C=Ci, tol=tolerance, cache_size=cachesize)
  else:
    wclf=SVC()


elif classificationtype == 'gaussianNB':
  from sklearn.naive_bayes import GaussianNB
  wclf=GaussianNB()
  
else:
  print 'no classifier chosen'
  sys.exit()

wclf.fit(traincounts,trainweight)


################ save the classifier ####################

import cPickle
with open('Classpredictor_'+classificationtype+os.path.splitext(output)[0]+'.pkl', 'wb') as fid:
    cPickle.dump(wclf, fid) 
  #np.savetxt('SVM-Parameter_'+outname, supvec , header=str(rcut)+' '+Benchmark+' '+kern+' probability: '+str(prob)+' tolerance: '+str(tolerance)+' C: '+str(Ci)) 

if predicton:
  testscore=wclf.score(testcounts,testweight)
  trainscore=wclf.score(traincounts,trainweight)
  print 'testaccuracy: ', testscore
  print 'trainingsaccuracy: ', trainscore
  if classificationtype == 'gaussianNB':
    prediction = -wclf.predict_proba(testcounts)[:,-1]
    trainprediction = -wclf.predict_proba(traincounts)[:,-1]
  else:
    prediction = wclf.decision_function(testcounts)
    trainprediction = wclf.decision_function(traincounts)
  np.savetxt('Classifier-Predictions_testset_'+classificationtype+output+'.dat',np.array(zip(prediction,testweight)),fmt="%6s")
  np.savetxt('Classifier-Predictions_training_'+classificationtype+output+'.dat',np.array(zip(trainprediction,trainweight)),fmt="%6s")
