import sys,os
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from matplotlib.widgets import Slider
from matplotlib import cm

if '--plotparams' in sys.argv:
  parname = sys.argv[sys.argv.index('--plotparams')+1]
  parameter = np.genfromtxt(parname,skip_header=1)
  atyps = len(parameter[0])
  if '--attractparameter' in sys.argv:
    power1 = float(sys.argv[sys.argv.index('--attractparameter')+1])
    power2 = float(sys.argv[sys.argv.index('--attractparameter')+2])
    Emin = ((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*parameter[atyps:2*atyps]*parameter[2*atyps:]
    rmin = ((power1/power2)**(1.0/(power1-power2)))*parameter[:atyps]
    parameter = np.append(rmin,Emin, axis = 0)
  
  if '--namelist' in sys.argv:
    names = np.genfromtxt(sys.argv[sys.argv.index('--namelist')+1], dtype = str)[:,-1]
  else:
    names = np.arange(1,len(parameter)+1,1).astype(str)
  
  if len(np.shape(parameter))>1:
    bins = len(parameter)/atyps
    sortpar = np.zeros((bins, (atyps*(atyps-1)/2)+atyps))
    sortnames = []
    ind = 0

    for n in range(atyps):
      for nn in range(n,atyps):
	for b in range(bins):
	  sortpar[b,ind] = parameter[b*atyps+n,nn]
	sortnames.append(names[n]+' / '+names[nn])
	ind += 1
    sortnames = np.array(sortnames, dtype = str)
    for b in range(bins):
      print 'bins ', b
      sort = np.argsort(sortpar[b])
      sortpar[b] = sortpar[b][sort]
      sortname = sortnames[sort]
      print ' unfavorable contacts               |                  favorable contacts '
      for i in range(30):
	print "%32s" % sortname[-i-1], "%3.3f" % sortpar[b][-i-1],"%32s" % sortname[i], "%3.3f" % sortpar[b][i]
      
      if '--sortplot' in sys.argv:
	fig=plt.figure('best')
	ax=fig.add_subplot(111)
	ax.set_ylabel('$\mathrm{p}$',rotation=0,fontsize=15)
	ax.xaxis.set_ticks_position('top')
	ax.set_xticks(np.arange(0,30))
	ax.set_xticklabels(sortname[:30],rotation=30 , ha='left')
	p=ax.bar(np.arange(0,30),sortpar[b][:30],0.5)
	fig1=plt.figure('worst')
	ax1=fig1.add_subplot(111)
	ax1.set_ylabel('$\mathrm{p}$',rotation=0, fontsize=15)
	ax1.set_xticks(np.arange(0,30))
	ax1.set_xticklabels(sortname[-30:], rotation =-30,ha='left')
	p1=ax1.bar(np.arange(0,30),sortpar[b][-30:][-np.arange(1,31)],0.5)
	plt.show()
	
	
	
    xnames = names
    ynames = []
    for i in range(bins):  
      ynames += names
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    p=ax.imshow(parameter,cmap='spectral',interpolation='none')
    ax.set_yticks(np.arange(0, len(parameter), 1))
    ax.set_xticks(np.arange(0, len(parameter[0]), 1))
    ax.set_xticklabels(xnames, rotation = 90)
    ax.set_yticklabels(ynames)
    fig.colorbar(p)
    plt.show()
  
  else:
    parameter = parameter.reshape((1,len(parameter)))
    xnames = names
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    p=ax.imshow(parameter,cmap='spectral',interpolation='none')
    ax.set_xticks(np.arange(0, len(parameter[0]), 1))
    ax.set_xticklabels(xnames, rotation = 90)
    fig.colorbar(p)
    plt.show()


if '--crossparameter' in sys.argv:
  numcross = int(sys.argv[sys.argv.index('--crossparameter')+1])
  parameter = []
  aton = False
  
  if '--attractparameter' in sys.argv:
    aton = True
    power1 = float(sys.argv[sys.argv.index('--attractparameter')+1])
    power2 = float(sys.argv[sys.argv.index('--attractparameter')+2])
  
  erron = False
  if '--meanerror' in sys.argv:
    erron = True
  
  for i in range(numcross):
    pars = np.genfromtxt(sys.argv[sys.argv.index('--crossparameter')+2+i])
    atyps=len(pars[0])
    if aton:
      Emin = ((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*pars[atyps:2*atyps]*pars[2*atyps:]
      rmin = ((power1/power2)**(1.0/(power1-power2)))*pars[:atyps]
      sign = pars[2*atyps:]
      pars = np.append(rmin,np.append(Emin,sign,axis=0), axis = 0)
    parameter.append(pars)
  
  
  parameter = np.array(parameter)
  
  shapepar = len(np.shape(parameter[0]))
  if shapepar > 1:
    bins=int(np.shape(parameter[0])[0]/atyps)      
    mean=np.zeros((numcross,bins))
    std=np.zeros((numcross,bins))
    for i in range(bins):
      for j in range(numcross):
	std[j,i]=np.std(parameter[j,i*atyps:(i+1)*atyps])
	if aton == False:
	  mean[j,i]=np.mean(parameter[j,i*atyps:(i+1)*atyps])
    
    if aton:
      std[:,2]= np.ones(numcross)
    
    corr=np.zeros(np.shape(parameter[0]))
    err = np.zeros(np.shape(parameter[0]))
    for b in range(bins):
      for n in range(atyps):
	for nn in range(atyps):
	  for i in range(numcross):
	    for j in range(i, numcross): 
	      corr[b*atyps+n,nn]+=(parameter[i,b*atyps+n,nn]-mean[i,b])*(parameter[j,b*atyps+n,nn]-mean[j,b])/(std[j,b]*std[i,b])
	      err[b*atyps+n,nn]+=((parameter[i,b*atyps+n,nn]-mean[i,b])-(parameter[j,b*atyps+n,nn]-mean[j,b]))**2/(std[j,b]*std[i,b])
      corr = corr/float(numcross)
      fig=plt.figure('Correlation')
      ax=fig.add_subplot(111)
      p=ax.imshow(corr[b*atyps:(b+1)*atyps],cmap='spectral',interpolation='none')
      fig.colorbar(p)
      
      if erron:
	err = err/float(numcross)
	fig1=plt.figure('Error')
	ax1=fig1.add_subplot(111)
	p1=ax1.imshow(err[b*atyps:(b+1)*atyps],cmap='spectral',interpolation='none')
	fig1.colorbar(p1)
      
      plt.show() 
   
  else:
    atyps=int(np.shape(parameter[0])[-1])    
    mean=np.zeros(numcross)
    std=np.zeros(numcross)
    for j in range(numcross):
	std[j]=np.std(parameter[j])
	mean[j]=np.mean(parameter[j])
    
    corr=np.zeros((1,len(parameter[0])))
    for n in range(atyps):
      for i in range(numcross):
	for j in range(i, numcross): 
	  corr[0,n]+=(parameter[i,n]-mean[i])*(parameter[j,n]-mean[j])/(std[j]*std[i])
    
    corr = corr/float(numcross)
    fig=plt.figure('Correlation')
    ax=fig.add_subplot(111)
    p=ax.imshow(corr,cmap='spectral',interpolation='none')
    fig.colorbar(p)
    plt.show()
    

if '--parameterdifference' in sys.argv:
  parameter1=np.genfromtxt(sys.argv[sys.argv.index('--parameterdifference')+1],skip_header=1)
  parameter2=np.genfromtxt(sys.argv[sys.argv.index('--parameterdifference')+2],skip_header=1)
  atyps=int(np.shape(parameter1)[-1])
  aton = False
  if '--attractparameter' in sys.argv:
    aton = True
    power1 = float(sys.argv[sys.argv.index('--attractparameter')+1])
    power2 = float(sys.argv[sys.argv.index('--attractparameter')+2])
    
    Emin = ((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*parameter1[atyps:2*atyps]*parameter1[2*atyps:]
    rmin = ((power1/power2)**(1.0/(power1-power2)))*parameter1[:atyps]
    sign = parameter1[2*atyps:]
    parameter1 = np.append(rmin,np.append(Emin,sign,axis=0), axis = 0)
  
    Emin = ((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*parameter2[atyps:2*atyps]*parameter2[2*atyps:]
    rmin = ((power1/power2)**(1.0/(power1-power2)))*parameter2[:atyps]
    sign = parameter2[2*atyps:]
    parameter2 = np.append(rmin,np.append(Emin,sign,axis=0), axis = 0)
  
  if np.shape(parameter1)==np.shape(parameter2):
    bins=int(np.shape(parameter1)[0]/atyps)
    shapepar = len(np.shape(parameter1))
    
    if aton:
      bins = 2
    
    if shapepar>1:
      mean=np.zeros((2,bins))
      std=np.zeros((2,bins))
      for i in range(bins):
	print 'for bin: ', i
	std[0,i]=np.std(parameter1[i*atyps:(i+1)*atyps])
	std[1,i]=np.std(parameter2[i*atyps:(i+1)*atyps])
	mean[0,i]=np.mean(parameter1[i*atyps:(i+1)*atyps])
	mean[1,i]=np.mean(parameter2[i*atyps:(i+1)*atyps])
	print 'max1', np.amax(parameter1[i*atyps:(i+1)*atyps])
	print 'min1', np.amin(parameter1[i*atyps:(i+1)*atyps])
	print 'max2', np.amax(parameter2[i*atyps:(i+1)*atyps])
	print 'min2', np.amin(parameter2[i*atyps:(i+1)*atyps])
	print 'mean1 ',mean[0,i]
	print 'mean2 ',mean[1,i]
	print 'std1 ',std[0,i]
	print 'std2 ',std[1,i]
	parameter1[i*atyps:(i+1)*atyps]=parameter1[i*atyps:(i+1)*atyps]/std[0,i]
	parameter2[i*atyps:(i+1)*atyps]=parameter2[i*atyps:(i+1)*atyps]/std[1,i]
    else:
      parameter1 = parameter1.reshape((1,atyps))
      parameter2 = parameter2.reshape((1,atyps))
      std1=np.std(parameter1)
      std2=np.std(parameter2)
      mean1=np.mean(parameter1)
      mean2=np.mean(parameter2)
      parameter1=parameter1/std1
      parameter2=parameter2/std2
      
    
    out=np.absolute(parameter1-parameter2)
    out2=np.ones(np.shape(out))-(np.absolute(parameter1-parameter2))/(np.absolute(parameter1)+np.absolute(parameter2))
    out3=parameter1*parameter2
    
    if shapepar > 1:
      for i in range(bins):
	parameter1[:(i+1)*atyps]=parameter1[:(i+1)*atyps]-mean[0,i]/std[0,i]
	parameter2[:(i+1)*atyps]=parameter2[:(i+1)*atyps]-mean[1,i]/std[1,i]
    else:
      parameter1=parameter1-mean1/std1
      parameter2=parameter2-mean2/std2
    
    out5=parameter1*parameter2
    corr=np.mean(out5)
    print 'corr', corr
    rmsd=np.mean(out**2)
    print 'rmsd ', rmsd
    
    fig=plt.figure('scaled difference')
    ax=fig.add_subplot(111)
    p=ax.imshow(out,cmap='spectral',interpolation='none')
    fig.colorbar(p)
    fig2=plt.figure('scaled 1-difference')
    ax2=fig2.add_subplot(111)
    p2=ax2.imshow(out2,cmap='spectral',interpolation='none')
    fig2.colorbar(p2)
    fig3=plt.figure('scaled product')
    ax3=fig3.add_subplot(111)
    p3=ax3.imshow(out3,cmap='spectral',interpolation='none')
    fig3.colorbar(p3)
    fig5=plt.figure('mean shifted product, corr')
    ax5=fig5.add_subplot(111)
    p5=ax5.imshow(out5,cmap='spectral',interpolation='none')
    fig5.colorbar(p5)
    plt.show()
  else:
    print 'parameter must have same shape!'
    sys.exit()
    

################## score comparison ####################




if '--scores' in sys.argv:
  numscores=int(sys.argv[sys.argv.index('--scores')+1])
  scorenames=[]
  for n in range(numscores):
    scorenames.append(sys.argv[sys.argv.index('--scores')+2+n])
    
  legendnames=['Real Ranking']
  if '--names' in sys.argv:
    nind = sys.argv.index('--names')
    for n, scorename in enumerate(scorenames):
      legendnames.append(sys.argv[nind+1+n])
      if scorename.split('.')[-1]=='dat':
	if '--vdwpotential' in sys.argv:
	  numscores+=1
	  legendnames.append('Vdw_'+sys.argv[nind+1+n])
	if '--electrostatics' in sys.argv:  
	  numscores+=1
	  legendnames.append('Elec_'+sys.argv[nind+1+n])
  else:  
    for scorename in scorenames:
      legendnames.append(scorename[:50])
      if scorename.split('.')[-1]=='dat':
	if '--vdwpotential' in sys.argv:
	  numscores+=1
	  legendnames.append('Vdw_'+scorename[:46])
	if '--electrostatics' in sys.argv:  
	  numscores+=1
	  legendnames.append('Elec_'+scorename[:45])
    
  dirlist = [x for x in os.listdir('.') if os.path.isdir(x)]

  benchon = False
  classon = False
  if '--benchmark' in sys.argv:
    setname=['trainingset', 'testset', 'None']
    benchon = True
    benchmark=np.genfromtxt(sys.argv[sys.argv.index('--benchmark')+1],dtype=str)
    sets = float(sys.argv[sys.argv.index('--benchmark')+2])
    if sets == 1:
      trainset = benchmark
    else:
      trainset=benchmark[:int(len(benchmark)*(sets-1.)/sets)]
    trainlen=len(trainset)
    for train in trainset:
      dirlist.remove(train)
    dirlist=np.append(trainset,dirlist)
    dirlist = np.array(dirlist)
    testlen=len(dirlist)-trainlen
    otherlen = 0
  elif '--classification' in sys.argv:
    classon = True
    classtype = sys.argv[sys.argv.index('--classification')+1]
    classfile = sys.argv[sys.argv.index('--classification')+2]
    cobj = open(classfile, 'r')
    com1=[]
    com2=[]
    com3=[]
    for cline in cobj:
      cline = cline.strip()
      cline = cline.split()
      com1.append(cline[0])
      com2.append(cline[1])
      com3.append(cline[2])
    cobj.close()
    comclass=[com1,com2,com3]
    if classtype == 'difficulty':
      setname=['rigid', 'medium', 'hard']
      classes = np.zeros(len(dirlist))
      dictio = { "rigid":0, "medium":1, "hard":2}
      for i, folder in enumerate(dirlist):
	classes[i] = dictio[comclass[1][comclass[0].index(folder)]]
    elif classtype == 'proteintype':
      setname=['enzyme', 'antibody', 'other']
      classes = np.zeros(len(dirlist))
      dictio = { "enzyme":0, "antibody":1, "other":2}
      for i, folder in enumerate(dirlist):
	classes[i] = dictio[comclass[2][comclass[0].index(folder)]]
    dirlist = np.array(dirlist)
    dirlist = dirlist[np.argsort(classes)]
    trainlen = len(ma.compressed(ma.masked_not_equal(classes,0)))
    testlen = len(ma.compressed(ma.masked_not_equal(classes,1)))
    otherlen = len(ma.compressed(ma.masked_not_equal(classes,2)))
  else:
    setname=['whole benchmark', 'None', 'None']
    trainlen=len(dirlist) 
    testlen = 0
    otherlen = 0
  
  print setname[0], trainlen
  print setname[1], testlen
  print setname[2], otherlen


  nstruc=100000 
  
  evaluations=[]
  if '--evaluate' in sys.argv:
    numscores+=1
    yvalues=sys.argv[sys.argv.index('--evaluate')+1]
    ytype=yvalues.split('.')[-1]
    ycut=float(sys.argv[sys.argv.index('--evaluate')+2])
    print 'read in yvalues ...'
    evals=[]
    yevals=[]
    sortevals=[]
    for folder in dirlist:
      yval=np.genfromtxt(folder+'/'+yvalues)
      if len(np.shape(yval))!=1:
	yval=yval[:,-1]
      if ytype[:8]=='capstars' or ytype[:8]=='fnat':
	probs=np.array(yval>=ycut,dtype=int)
	sort=np.argsort(-yval)
      if ytype[:5]=='irmsd' or ytype[:5]=='lrmsd':
	probs=np.array(yval<=ycut,dtype=int)
	sort=np.argsort(yval)
      evals.append(probs)
      yevals.append(yval)
      sortevals.append(sort)
    evaluations.append(sortevals)
    
  else:
    print 'please insert values to evaluate score on!'
    sys.exit()  

  for scorename in scorenames:
    elon=True
    vdwon=True
    if scorename.split('.')[-1]=='dat':
      if '--electrostatics' in sys.argv:
	elon=False
	if os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_el.npz'):
	  elon=True
      if '--vdwpotential' in sys.argv:
	vdwon=False
	if os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_vdw.npz'):
	  vdwon=True	
    if vdwon and elon and os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'.npz'):
      print 'load Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'.npz ...'
      inside=np.load('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'.npz')
      folderorder=inside['dirlist'].tolist()
      evaluation=inside['evaluation']
      folderindex=[]
      for folder in dirlist:
	folderindex.append(folderorder.index(folder))
      evaluations.append(evaluation[folderindex])
      if '--vdwpotential' in sys.argv:
	if os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_vdw.npz'):
	  inside=np.load('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_vdw.npz')
	  evaluation=inside['vdwevaluation']
	  folderorder=inside['dirlist'].tolist()
	  folderindex=[]
	  for folder in dirlist:
	    folderindex.append(folderorder.index(folder))
	  evaluations.append(evaluation[folderindex])
      if '--electrostatics' in sys.argv:
	if os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_el.npz'):
	  inside=np.load('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_el.npz')
	  evaluation=inside['elevaluation']
	  folderorder=inside['dirlist'].tolist()
	  folderindex=[]
	  for folder in dirlist:
	    folderindex.append(folderorder.index(folder))
	  evaluations.append(evaluation[folderindex])
      print '...load done'
    else:
      print '...read in ', scorename
      evaluation=[]
      vdwevaluation=[]
      elevaluation=[]
      elon=False
      vdwon=False
      for f,folder in enumerate(dirlist):
	if scorename.split('.')[-1]=='rescore':
	  scoring=np.genfromtxt(folder+'/'+scorename)
	elif scorename.split('.')[-1]=='dat':
	  scoring=[]
	  elec=np.empty(0)
	  vdws=np.empty(0)
	  gobj=open(folder+'/'+scorename,'r')
	  glines=gobj.readlines()
	  for g,gline in enumerate(glines):
	    if gline[:10]=='## Energy:':
	      gline=gline.strip()
	      gline=gline.split()
	      scoring.append(float(gline[2]))
	      if '--electrostatics' in sys.argv:
		ele=glines[g+1].strip()
		ele=ele.split()
		elec=np.append(elec,[float(ele[2])])
	      if '--vdwpotential' in sys.argv:
		v=glines[g+1].strip()
		v=v.split()
		if v[1][0]=='*':
		  vdws=np.append(vdws,[10000000])
		else:
		  vdws=np.append(vdws,[float(v[1])])
	  if '--vdwpotential' in sys.argv:
	    vdwon=True
	    vdws=np.array(vdws)
	    vsort=np.argsort(vdws)
	    vdwevaluation.append(vsort)
	  if '--electrostatics' in sys.argv:
	    elon=True
	    elec=np.array(elec)
	    esort=np.argsort(elec)
	    elevaluation.append(esort)
	else:
	  print 'it must be a file .rescore or .dat'
	  sys.exit()
	scoring=np.array(scoring)
	sort=np.argsort(scoring)
	evaluation.append(sort)
      evaluations.append(evaluation)
      
      np.savez('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'.npz',dirlist=dirlist,evaluation=evaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
      if vdwon:
	evaluations.append(vdwevaluation)
	np.savez('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_vdw.npz',dirlist=dirlist,vdwevaluation=vdwevaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
      if elon:
	evaluations.append(elevaluation)
	np.savez('Argsort_'+os.path.splitext(scorename)[0]+'-'+ytype+'-'+str(ycut)+'_el.npz',dirlist=dirlist,elevaluation=elevaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
      
  evaluations=np.array(evaluations)
  
  
  if '--Rankcorrstructures' in sys.argv:

    if '--goodstructures' in sys.argv:
      print 'only goodstructures'
      corr=np.ones((len(dirlist),numscores,numscores))
      for count in range(len(dirlist)):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  mean=np.mean(evaluations[0][count])
	  std=np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      corr[count,n,nn]=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)
	      corr[count,nn,n]=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)
	else:
	  corr[count]=np.eye(numscores)
    else:
      corr=np.ones((len(dirlist),numscores,numscores))
      for count in range(len(dirlist)):
	mean=np.mean(evaluations[0][count])
	std=np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    corr[count,n,nn]=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)
	    corr[count,nn,n]=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)

    fig = plt.figure('Rankcorrelation_complexes')
    ax = fig.add_subplot(111)
    frame = 0
    ax.set_title(dirlist[frame])
    ax.set_yticks(np.arange(0,numscores,1))
    ax.set_yticklabels(legendnames)
    ax.set_xticks(np.arange(0,numscores,1))
    ax.set_xticklabels(legendnames, rotation=60)
    ln = ax.imshow(corr[frame],cmap='spectral',interpolation='none')
    fig.colorbar(ln)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Frame', 0, len(dirlist)-1, valinit=0,valfmt='%d')
    # call back function
    def update(val):
      frame = np.floor(sframe.val)
      ln.set_data(corr[frame])
      ax.set_title(dirlist[frame])
      plt.draw() 
    sframe.on_changed(update)
    plt.show()


  if '--Rankerrorstructures' in sys.argv:

    if '--goodstructures' in sys.argv:
      print 'only goodstructures'
      err=np.zeros((len(dirlist),numscores,numscores))
      for count in range(len(dirlist)):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  std=np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      err[count,n,nn]=np.sum((ma.compressed(ma.array(evaluations[n][count]-evaluations[nn][count], mask = maske)))**2/std**2)
	      err[count,nn,n]=np.sum((ma.compressed(ma.array(evaluations[nn][count]-evaluations[n][count], mask = maske)))**2/std**2)

	else:
	  err[count]=np.ones((numscores,numscores))-np.eye(numscores)
    else:
      err=np.zeros((len(dirlist),numscores,numscores))
      for count in range(len(dirlist)):
	std=np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    err[count,n,nn]=np.sum((evaluations[n][count]-evaluations[nn][count])**2/std**2)
	    err[count,nn,n]=np.sum((evaluations[nn][count]-evaluations[n][count])**2/std**2)

    fig = plt.figure('Rankerror_complexes')
    ax = fig.add_subplot(111)
    frame = 0
    ax.set_title(dirlist[frame])
    ax.set_yticks(np.arange(0,numscores,1))
    ax.set_yticklabels(legendnames)
    ax.set_xticks(np.arange(0,numscores,1))
    ax.set_xticklabels(legendnames, rotation=60)
    ln = ax.imshow(err[frame],cmap='spectral',interpolation='none')
    fig.colorbar(ln)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Frame', 0, len(dirlist)-1, valinit=0,valfmt='%d')
    # call back function
    def update(val):
      frame = np.floor(sframe.val)
      ln.set_data(err[frame])
      ax.set_title(dirlist[frame])
      plt.draw() 
    sframe.on_changed(update)
    plt.show()
    
  
  if '--Rankerror' in sys.argv:
    if '--goodstructures' in sys.argv:
      err = np.zeros((3, numscores, numscores))
      for count in range(trainlen):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  std = np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      err[0,n,nn] += np.sum((ma.compressed(ma.array(evaluations[n][count]-evaluations[nn][count], mask = maske)))**2/std**2)/float(trainlen)
	      err[0,nn,n] += np.sum((ma.compressed(ma.array(evaluations[nn][count]-evaluations[n][count], mask = maske)))**2/std**2)/float(trainlen)
      for count in range(trainlen, trainlen +testlen):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  std = np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      err[1,n,nn] += np.mean((ma.compressed(ma.array(evaluations[n][count]-evaluations[nn][count], mask = maske)))**2/std**2)/float(testlen)
	      err[1,nn,n] += np.mean((ma.compressed(ma.array(evaluations[nn][count]-evaluations[n][count], mask = maske)))**2/std**2)/float(testlen)  
      for count in range(trainlen+testlen, len(dirlist)):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  std = np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      err[2,n,nn] += np.mean((ma.compressed(ma.array(evaluations[n][count]-evaluations[nn][count], mask = maske)))**2/std**2)/float(otherlen)
	      err[2,nn,n] += np.mean((ma.compressed(ma.array(evaluations[nn][count]-evaluations[n][count], mask = maske)))**2/std**2)/float(otherlen)
    else:
      err = np.zeros((3, numscores, numscores))
      for count in range(trainlen):
	std = np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    err[0,n,nn] += np.mean((evaluations[n][count]-evaluations[nn][count])**2/std**2)/float(trainlen)
	    err[0,nn,n] += np.mean((evaluations[nn][count]-evaluations[n][count])**2/std**2)/float(trainlen)
      for count in range(trainlen, trainlen +testlen):
	std = np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    err[1,n,nn] += np.mean((evaluations[n][count]-evaluations[nn][count])**2/std**2)/float(testlen)
	    err[1,nn,n] += np.mean((evaluations[nn][count]-evaluations[n][count])**2/std**2)/float(testlen)  
      for count in range(trainlen+testlen, len(dirlist)):
	std = np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    err[2,n,nn] += np.mean((evaluations[n][count]-evaluations[nn][count])**2/std**2)/float(otherlen)
	    err[2,nn,n] += np.mean((evaluations[nn][count]-evaluations[n][count])**2/std**2)/float(otherlen)  
  
    print err
    fig=plt.figure('Rankerror_'+setname[0]+'_'+str(trainlen)+'_complexes')
    ax=fig.add_subplot(111)
    ax.set_yticks(np.arange(0,numscores,1))
    ax.set_yticklabels(legendnames)
    ax.set_xticks(np.arange(0,numscores,1))
    ax.set_xticklabels(legendnames, rotation=90)
    p=ax.imshow(err[0],cmap='spectral',interpolation='none')
    fig.colorbar(p)
    if benchon or classon:
      fig2=plt.figure('Rankerror_'+setname[1]+'_'+str(testlen)+'_complexes')
      ax2=fig2.add_subplot(111)
      ax2.set_yticks(np.arange(0,numscores,1))
      ax2.set_yticklabels(legendnames)
      ax2.set_xticks(np.arange(0,numscores,1))
      ax2.set_xticklabels(legendnames, rotation=90)
      p2=ax2.imshow(err[1],cmap='spectral',interpolation='none')
      fig2.colorbar(p)
    if classon:
      fig3=plt.figure('Rankerror_'+setname[2]+'_'+str(otherlen)+'_complexes')
      ax3=fig3.add_subplot(111)
      ax3.set_yticks(np.arange(0,numscores,1))
      ax3.set_yticklabels(legendnames)
      ax3.set_xticks(np.arange(0,numscores,1))
      ax3.set_xticklabels(legendnames, rotation=90)
      p3=ax3.imshow(err[2],cmap='spectral',interpolation='none')
      fig3.colorbar(p)

    plt.show()
    
    
  if '--Rankcorr' in sys.argv:
    if '--goodstructures' in sys.argv:
      print 'only goodstructures'
      corr=np.zeros((3,numscores,numscores))+np.eye(numscores)
      for count in range(trainlen):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  mean=np.mean(evaluations[0][count])
	  std=np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      corr[0,n,nn]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/trainlen
	      corr[0,nn,n]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/trainlen

      for count in range(trainlen, trainlen+testlen):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  mean=np.mean(evaluations[0][count])
	  std=np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      corr[1,n,nn]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/testlen
	      corr[1,nn,n]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/testlen
      for count in range(trainlen+testlen, len(dirlist)):
	maske=ma.getmask(ma.masked_equal(evals[count],0))
	if len(ma.nonzero(evals[count])[0])!=0:
	  mean=np.mean(evaluations[0][count])
	  std=np.std(evaluations[0][count])
	  for n in range(numscores):
	    for nn in range(n+1,numscores):
	      corr[2,n,nn]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/testlen
	      corr[2,nn,n]+=np.mean((ma.compressed(ma.array(evaluations[n][count],mask=maske))-mean)/std*(ma.compressed(ma.array(evaluations[nn][count],mask=maske))-mean)/std)/testlen

    else:
      corr=np.zeros((3,numscores,numscores))+np.eye(numscores)
      for count in range(trainlen):
	mean=np.mean(evaluations[0][count])
	std=np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    corr[0,n,nn]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/trainlen
	    corr[0,nn,n]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/trainlen
      for count in range(trainlen, trainlen+testlen):
	mean=np.mean(evaluations[0][count])
	std=np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    corr[1,n,nn]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/testlen
	    corr[1,nn,n]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/testlen
      for count in range(trainlen+testlen, len(dirlist)):
	mean=np.mean(evaluations[0][count])
	std=np.std(evaluations[0][count])
	for n in range(numscores):
	  for nn in range(n+1,numscores):
	    corr[2,n,nn]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/testlen
	    corr[2,nn,n]+=np.mean((evaluations[n][count]-mean)/std*(evaluations[nn][count]-mean)/std)/testlen
        
    print corr
    fig=plt.figure('Rankcorrelation_'+setname[0]+'_'+str(trainlen)+'_complexes')
    ax=fig.add_subplot(111)
    ax.set_yticks(np.arange(0,numscores,1))
    ax.set_yticklabels(legendnames)
    ax.set_xticks(np.arange(0,numscores,1))
    ax.set_xticklabels(legendnames, rotation=60)
    p=ax.imshow(corr[0],cmap='spectral',interpolation='none')
    fig.colorbar(p)
    if benchon or classon:
      fig2=plt.figure('Rankcorrelation_'+setname[1]+'_'+str(testlen)+'_complexes')
      ax2=fig2.add_subplot(111)
      ax2.set_yticks(np.arange(0,numscores,1))
      ax2.set_yticklabels(legendnames)
      ax2.set_xticks(np.arange(0,numscores,1))
      ax2.set_xticklabels(legendnames, rotation=60)
      p2=ax2.imshow(corr[1],cmap='spectral',interpolation='none')
      fig2.colorbar(p)
    if classon:
      fig3=plt.figure('Rankcorrelation_'+setname[2]+'_'+str(otherlen)+'_complexes')
      ax3=fig3.add_subplot(111)
      ax3.set_yticks(np.arange(0,numscores,1))
      ax3.set_yticklabels(legendnames)
      ax3.set_xticks(np.arange(0,numscores,1))
      ax3.set_xticklabels(legendnames, rotation=60)
      p3=ax3.imshow(corr[2],cmap='spectral',interpolation='none')
      fig3.colorbar(p)

    plt.show()
  
  if '--topcorr' in sys.argv:
    top=float(sys.argv[sys.argv.index('--topcorr')+1])
    toptype = sys.argv[sys.argv.index('--topcorr')+2]
    symdif=np.zeros((3,numscores,numscores))
    union=np.zeros((3,numscores,numscores))
    complement=np.zeros((3,numscores,numscores))
    symon=False
    unionon=False
    comon=False
    if '--symdifference' in sys.argv:
      symon=True
    if '--union' in sys.argv:
      unionon=True
    if '--complement' in sys.argv:
      comon=True
    
    if unionon:
      if toptype == 'all':
	for z in range(3):
	  union[z]=np.eye(np.shape(union)[1])
      if toptype == 'bad':
	for z in range(3):
	  union[z]=np.eye(np.shape(union)[1])
      if toptype == 'good':
	for n in range(numscores):
	  for count in range(trainlen):
	    ind=int(top*len(evals[count]))
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)	 
	    sortstruc1=struc[evaluations[n][count]]
	    struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	    union[0,n,n]+=len(struc1)/float(maxstruc*trainlen)
	  for count in range(trainlen,trainlen+testlen):
	    ind=int(top*len(evals[count]))
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)	 
	    sortstruc1=struc[evaluations[n][count]]
	    struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	    union[1,n,n]+=len(struc1)/float(maxstruc*testlen)
	  for count in range(trainlen+testlen, len(dirlist)):
	    ind=int(top*len(evals[count]))
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)	 
	    sortstruc1=struc[evaluations[n][count]]
	    struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	    union[2,n,n]+=len(struc1)/float(maxstruc*otherlen)
	    
	    
    for n in range(numscores):
      for nn in range(n+1,numscores):
	for count in range(trainlen):
	  ind=int(top*len(evals[count]))
	  if toptype == 'all':
	    maxstruc = ind
	    struc = np.arange(1,len(evals[count])+1)
	  elif toptype == 'good':
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)	    
	  elif toptype == 'bad':
	    compevals = np.absolute(evals[count]-1.)
	    maxstruc = ind
	    struc=compevals*np.arange(1,len(compevals)+1)
	  else:
	    print 'toptype not understood'
	    sys.exit()

	  sortstruc1=struc[evaluations[n][count]]
	  sortstruc2=struc[evaluations[nn][count]]
	  struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	  struc2=ma.compressed(ma.masked_equal(sortstruc2[:ind],0))
	  if symon:
	    sett=np.setxor1d(struc1,struc2)
	    symdif[0,n,nn]+=len(sett)/float(maxstruc*trainlen)
	    symdif[0,nn,n]+=len(sett)/float(maxstruc*trainlen)
	  if unionon:
	    sett=np.union1d(struc1,struc2)
	    union[0,n,nn]+=len(sett)/float(maxstruc*trainlen)
	    union[0,nn,n]+=len(sett)/float(maxstruc*trainlen)
	  if comon:
	    sett=np.setdiff1d(struc1,struc2)
	    sett2=np.setdiff1d(struc2,struc1)
	    complement[0,n,nn]+=len(sett)/float(maxstruc*trainlen)
	    complement[0,nn,n]+=len(sett2)/float(maxstruc*trainlen)
	for count in range(trainlen, trainlen+testlen):
	  ind=int(top*len(evals[count]))
	  if toptype == 'all':
	    maxstruc = ind
	    struc = np.arange(1,len(evals[count])+1)
	  elif toptype == 'good':
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)
	  elif toptype == 'bad':
	    compevals = np.absolute(evals[count]-1.)
	    maxstruc=ind
	    struc=compevals*np.arange(1,len(compevals)+1)
	  else:
	    print 'toptype not understood'
	    sys.exit()

	  sortstruc1=struc[evaluations[n][count]]
	  sortstruc2=struc[evaluations[nn][count]]

	  struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	  struc2=ma.compressed(ma.masked_equal(sortstruc2[:ind],0))
	  if symon:
	    sett=np.setxor1d(struc1,struc2)
	    symdif[1,n,nn]+=len(sett)/float(maxstruc*testlen)
	    symdif[1,nn,n]+=len(sett)/float(maxstruc*testlen)
	  if unionon:
	    sett=np.union1d(struc1,struc2)
	    union[1,n,nn]+=len(sett)/float(maxstruc*testlen)
	    union[1,nn,n]+=len(sett)/float(maxstruc*testlen)
	  if comon:
	    sett=np.setdiff1d(struc1,struc2)
	    sett2=np.setdiff1d(struc2,struc1)
	    complement[1,n,nn]+=len(sett)/float(maxstruc*testlen)
	    complement[1,nn,n]+=len(sett2)/float(maxstruc*testlen)
	for count in range(trainlen + testlen, len(dirlist)):
	  ind=int(top*len(evals[count]))
	  if toptype == 'all':
	    maxstruc = ind
	    struc = np.arange(1,len(evals[count])+1)
	  elif toptype == 'good':
	    maxstruc=max([np.sum(evals[count]),1])
	    struc=evals[count]*np.arange(1,len(evals[count])+1)
	  elif toptype == 'bad':
	    compevals = np.absolute(evals[count]-1.)
	    maxstruc= ind
	    struc=compevals*np.arange(1,len(compevals)+1)
	  else:
	    print 'toptype not understood'
	    sys.exit()

	  sortstruc1=struc[evaluations[n][count]]
	  sortstruc2=struc[evaluations[nn][count]]

	  struc1=ma.compressed(ma.masked_equal(sortstruc1[:ind],0))
	  struc2=ma.compressed(ma.masked_equal(sortstruc2[:ind],0))
	  if symon:
	    sett=np.setxor1d(struc1,struc2)
	    symdif[2,n,nn]+=len(sett)/float(maxstruc*otherlen)
	    symdif[2,nn,n]+=len(sett)/float(maxstruc*otherlen)
	  if unionon:
	    sett=np.union1d(struc1,struc2)
	    union[2,n,nn]+=len(sett)/float(maxstruc*otherlen)
	    union[2,nn,n]+=len(sett)/float(maxstruc*otherlen)
	  if comon:
	    sett=np.setdiff1d(struc1,struc2)
	    sett2=np.setdiff1d(struc2,struc1)
	    complement[2,n,nn]+=len(sett)/float(maxstruc*otherlen)
	    complement[2,nn,n]+=len(sett2)/float(maxstruc*otherlen)
    if symon:
      fig=plt.figure('symdifference_'+setname[0]+'_'+str(trainlen)+'_complexes_top_'+str(top))
      ax=fig.add_subplot(111)
      ax.set_yticks(np.arange(0,numscores,1))
      ax.set_yticklabels(legendnames[1:])
      ax.set_xticks(np.arange(0,numscores,1))
      ax.set_xticklabels(legendnames[1:], rotation=90)
      p=ax.imshow(symdif[0,1:,1:],cmap='spectral',interpolation='none')
      fig.colorbar(p)
      if benchon or classon:
	fig2=plt.figure('symdifference_'+setname[1]+'_'+str(testlen)+'_complexes_top_'+str(top))
	ax2=fig2.add_subplot(111)
	ax2.set_yticks(np.arange(0,numscores,1))
	ax2.set_yticklabels(legendnames[1:])
	ax2.set_xticks(np.arange(0,numscores,1))
	ax2.set_xticklabels(legendnames[1:], rotation=90)
	p2=ax2.imshow(symdif[1,1:,1:],cmap='spectral',interpolation='none')
	fig2.colorbar(p2)
      if classon:
	fig7=plt.figure('symdifference_'+setname[2]+'_'+str(otherlen)+'_complexes_top_'+str(top))
	ax7=fig7.add_subplot(111)
	ax7.set_yticks(np.arange(0,numscores,1))
	ax7.set_yticklabels(legendnames[1:])
	ax7.set_xticks(np.arange(0,numscores,1))
	ax7.set_xticklabels(legendnames[1:], rotation=90)
	p7=ax7.imshow(symdif[2,1:,1:],cmap='spectral',interpolation='none')
	fig7.colorbar(p7)
    if unionon:
      fig3=plt.figure('union_'+setname[0]+'_'+str(trainlen)+'_complexes_top_'+str(top))
      ax3=fig3.add_subplot(111)
      ax3.set_yticks(np.arange(0,numscores,1))
      ax3.set_yticklabels(legendnames[1:])
      ax3.set_xticks(np.arange(0,numscores,1))
      ax3.set_xticklabels(legendnames[1:], rotation=90)
      p3=ax3.imshow(union[0,1:,1:],cmap='spectral',interpolation='none')
      fig3.colorbar(p3)
      if benchon or classon:
	fig4=plt.figure('union_'+setname[1]+'_'+str(testlen)+'_complexes_top_'+str(top))
	ax4=fig4.add_subplot(111)
	ax4.set_yticks(np.arange(0,numscores,1))
	ax4.set_yticklabels(legendnames[1:])
	ax4.set_xticks(np.arange(0,numscores,1))
	ax4.set_xticklabels(legendnames[1:], rotation=90)
	p4=ax4.imshow(union[1,1:,1:],cmap='spectral',interpolation='none')
	fig4.colorbar(p4)
      if classon:
	fig8=plt.figure('union_'+setname[2]+'_'+str(otherlen)+'_complexes_top_'+str(top))
	ax8=fig8.add_subplot(111)
	ax8.set_yticks(np.arange(0,numscores,1))
	ax8.set_yticklabels(legendnames[1:])
	ax8.set_xticks(np.arange(0,numscores,1))
	ax8.set_xticklabels(legendnames[1:], rotation=90)
	p8=ax8.imshow(union[2,1:,1:],cmap='spectral',interpolation='none')
	fig8.colorbar(p8)
    if comon:
      fig6=plt.figure('complement_'+setname[0]+'_'+str(trainlen)+'_complexes_top_'+str(top))
      ax6=fig6.add_subplot(111)
      ax6.set_yticks(np.arange(0,numscores,1))
      ax6.set_yticklabels(legendnames[1:])
      ax6.set_xticks(np.arange(0,numscores,1))
      ax6.set_xticklabels(legendnames[1:], rotation=90)
      p6=ax6.imshow(complement[0,1:,1:],cmap='spectral',interpolation='none')
      fig6.colorbar(p6)
      if benchon or classon:
	fig5=plt.figure('complement_'+setname[1]+'_'+str(testlen)+'_complexes_top_'+str(top))
	ax5=fig5.add_subplot(111)
	ax5.set_yticks(np.arange(0,numscores,1))
	ax5.set_yticklabels(legendnames[1:])
	ax5.set_xticks(np.arange(0,numscores,1))
	ax5.set_xticklabels(legendnames[1:],rotation=90)
	p5=ax5.imshow(complement[1,1:,1:],cmap='spectral',interpolation='none')
	fig5.colorbar(p5)
      if classon:
	fig9=plt.figure('complement_'+setname[2]+'_'+str(testlen)+'_complexes_top_'+str(top))
	ax9=fig9.add_subplot(111)
	ax9.set_yticks(np.arange(0,numscores,1))
	ax9.set_yticklabels(legendnames[1:])
	ax9.set_xticks(np.arange(0,numscores,1))
	ax9.set_xticklabels(legendnames[1:], rotation=90)
	p9=ax9.imshow(complement[2,1:,1:],cmap='spectral',interpolation='none')
	fig9.colorbar(p9)
    plt.show()
    
if '--topplot' in sys.argv:
  scorename = sys.argv[sys.argv.index('--topplot')+1]
  strucident = sys.argv[sys.argv.index('--topplot')+2]
  xplot = sys.argv[sys.argv.index('--topplot')+3]
  yplot = sys.argv[sys.argv.index('--topplot')+4]
  top = int(sys.argv[sys.argv.index('--topplot')+5])
  
  dirlist = [x for x in os.listdir('.') if os.path.isdir(x)]
   
  if os.path.isfile('Argsort_'+os.path.splitext(scorename)[0]+'.npy'):
    print 'load...'
    evaluation=np.load('Argsort_'+os.path.splitext(scorename)[0]+'.npy')
  else:
    print 'evaluate...'
    evaluation=[]
    for f,folder in enumerate(dirlist):
      if scorename.split('.')[-1]=='rescore':
	scoring=np.genfromtxt(folder+'/'+scorename)
      elif scorename.split('.')[-1]=='dat':
	scoring=[]
	gobj=open(folder+'/'+scorename,'r')
	glines=gobj.readlines()
	for g,gline in enumerate(glines):
	  if gline[:10]=='## Energy:':
	    gline=gline.strip()
	    gline=gline.split()
	    scoring.append(float(gline[2]))
      else:
	print 'it must be a file .rescore or .dat'
	sys.exit()
      scoring=np.array(scoring)
      sort=np.argsort(scoring)
      evaluation.append(sort)
    np.save('Argsort_'+os.path.splitext(scorename)[0]+'.npy',evaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
  if '--plot' in sys.argv:
    print 'read in yvalues ...'
    identevals=[]
    xevals=[]
    yevals=[]
    for f,folder in enumerate(dirlist):
      print folder
      ival=np.genfromtxt(folder+'/'+strucident)
      xval=np.genfromtxt(folder+'/'+xplot)
      yval=np.genfromtxt(folder+'/'+yplot)
      if len(np.shape(xval))!=1:
	  xval=xval[:,-1]
      if len(np.shape(yval))!=1:
	  yval=yval[:,-1]
      if len(np.shape(ival))!=1:
	  ival=ival[:,-1]
      identevals.append(ival[evaluation[f]])
      xevals.append(xval[evaluation[f]])
      yevals.append(yval[evaluation[f]])
    
    plotoutgood=np.empty((2,0))
    plotoutneg=np.empty((2,0))
    for f,folder in enumerate(dirlist):
      maskegood=ma.getmask(ma.masked_less(identevals[f][:top],1))
      maskeneg=ma.getmask(ma.masked_less(identevals[f][top:],1))
      if len(ma.compressed(ma.masked_less(identevals[f][:top],1)))!=0:
	plotoutgood=[np.append(plotoutgood[0],ma.compressed(ma.array(xevals[f][:top],mask=maskegood))),np.append(plotoutgood[1],ma.compressed(ma.array(yevals[f][:top],mask=maskegood)))]
	plotoutneg=[np.append(plotoutneg[0],ma.compressed(ma.array(xevals[f][top:],mask=maskeneg))),np.append(plotoutneg[1],ma.compressed(ma.array(yevals[f][top:],mask=maskeneg)))]  

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel('$\mathrm{fnat}$', fontsize=15)
    ax.set_ylabel('$\mathrm{lrmsd}$',fontsize=15,rotation=0)
    p=ax.plot(plotoutneg[0],plotoutneg[1],'bo',markersize=7)
    p2=ax.plot(plotoutgood[0],plotoutgood[1],'ro', markersize=7)
    plt.show()
    
  
  
  
