# comparescores.py
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy as np
import sys,os
from matplotlib.widgets import Slider
from matplotlib import cm

dirlist = [x for x in os.listdir('.') if os.path.isdir(x)]

if '--scores' in sys.argv:
  numscores=int(sys.argv[sys.argv.index('--scores')+1])
  scorenames=[]
  for n in range(numscores):
    scorenames.append(sys.argv[sys.argv.index('--scores')+2+n])
else:
  print 'please insert scores to evaluate!'
  sys.exit()
  
legendnames=[]  
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
      legendnames.append(scorename)
      if scorename.split('.')[-1]=='dat':
	if '--vdwpotential' in sys.argv:
	  numscores+=1
	  legendnames.append('Vdw_'+scorename)
	if '--electrostatics' in sys.argv:  
	  numscores+=1
	  legendnames.append('Elec_'+scorename)

benchon = False
classon = False
if '--benchmark' in sys.argv:
  benchon=True
  benchmark=np.genfromtxt(sys.argv[sys.argv.index('--benchmark')+1],dtype=str)
  sets = float(sys.argv[sys.argv.index('--benchmark')+2])
if '--classification' in sys.argv:
  classon = True
  classtype = sys.argv[sys.argv.index('--classification')+1]
  classfile = sys.argv[sys.argv.index('--classification')+2]


if benchon and classon== False:
  setname=['trainingset', 'testset', 'None']
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
  otherlen= 0
  classes = np.zeros(len(dirlist))
elif classon and benchon == False:
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
  classes = np.zeros(len(dirlist))
elif classon and benchon:
  setname=['trainingset', 'testset', 'None']
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
  otherlen= 0
  
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
    classname=['rigid', 'medium', 'hard']
    classes = np.zeros(len(dirlist))
    dictio = { "rigid":0, "medium":1, "hard":2}
    for i, folder in enumerate(dirlist):
      classes[i] = dictio[comclass[1][comclass[0].index(folder)]]
  elif classtype == 'proteintype':
    classname=['enzyme', 'antibody', 'other']
    classes = np.zeros(len(dirlist))
    dictio = { "enzyme":0, "antibody":1, "other":2}
    for i, folder in enumerate(dirlist):
      classes[i] = dictio[comclass[2][comclass[0].index(folder)]]
  dirlist = np.array(dirlist)
  classlen0 = len(ma.compressed(ma.masked_not_equal(classes,0)))
  classlen1 = len(ma.compressed(ma.masked_not_equal(classes,1)))
  classlen2 = len(ma.compressed(ma.masked_not_equal(classes,2)))
else:
  setname=['whole benchmark', 'None', 'None']
  trainlen=len(dirlist) 
  testlen = 0
  otherlen = 0
  
nstruc=100000

if '--evaluate' in sys.argv:
  yvalues=sys.argv[sys.argv.index('--evaluate')+1]
  ytype=yvalues.split('.')[-1]
  ycutval=sys.argv[sys.argv.index('--evaluate')+2].replace('.','')
  probline=''
  if ycutval=='probas':
    probavalues=np.genfromtxt(sys.argv[sys.argv.index('probas')+1])
    probline=int(sys.argv[sys.argv.index('probas')+2])
  check=0
  for n, scorename in enumerate(scorenames):
    if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_'+ytype+'_'+ycutval+str(probline)+'.npz'):
      check+=1
    if scorename.split('.')[-1]=='dat':
      if '--vdwpotential' in sys.argv:
	if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_vdw_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	  check+=1
      if '--electrostatics' in sys.argv:
	if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_el_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	  check+=1
	  
  if check!=numscores: 
    print 'read in yvalues ...'
    evals=[]
    setlength=[]
    for folder in dirlist:
      yval=np.genfromtxt(folder+'/'+yvalues)
      if len(np.shape(yval))!=1:
	yval=yval[:,-1]
      setlength.append(len(yval))
      if ycutval=='probas':
	print 'giving probabilities...'
	probavalues=np.genfromtxt(sys.argv[sys.argv.index('probas')+1])
	probline=int(sys.argv[sys.argv.index('probas')+2])
	probs=np.zeros(np.shape(yval))
	for i in range(len(probavalues)-1):
	    rmin=probavalues[i,0]
	    rmax=probavalues[i+1,0]
	    maske=ma.masked_outside(yval,rmin,rmax)
	    coor=ma.nonzero(maske)
	    probs[coor[0]]=probavalues[i,probline]
	evals.append(probs)
      else:
	ycut=float(sys.argv[sys.argv.index('--evaluate')+2])
	if ytype[:8]=='capstars' or ytype[:4]=='fnat':
	  probs=np.array(yval>=ycut,dtype=int)
	if ytype[:5]=='irmsd' or ytype[:5]=='lrmsd':
	  probs=np.array(yval<=ycut,dtype=int)
	evals.append(probs)

else:
  print 'please insert values to evaluate score on!'
  sys.exit()
  
print '...read in done'

evaluations=[]
for n, scorename in enumerate(scorenames):
  elon=True
  vdwon=True
  if scorename.split('.')[-1]=='dat':
    if '--electrostatics' in sys.argv:
      elon=False
      if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_el_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	elon=True
    if '--vdwpotential' in sys.argv:
      vdwon=False
      if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_vdw_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	vdwon=True	
  if vdwon and elon and os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_'+ytype+'_'+ycutval+str(probline)+'.npz'):
    print 'load Evaluation_'+os.path.splitext(scorename)[0]+'_'+ytype+'_'+ycutval+str(probline)+'.npz ...'
    inside=np.load('Evaluation_'+os.path.splitext(scorename)[0]+'_'+ytype+'_'+ycutval+str(probline)+'.npz')
    folderorder=inside['dirlist'].tolist()
    evaluation=inside['evaluation']
    folderindex=[]
    for folder in dirlist:
      folderindex.append(folderorder.index(folder))
    evaluations.append(evaluation[folderindex])
    if '--vdwpotential' in sys.argv:
      if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_vdw_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	inside=np.load('Evaluation_'+os.path.splitext(scorename)[0]+'_vdw_'+ytype+'_'+ycutval+str(probline)+'.npz')
	evaluation=inside['vdwevaluation']
	folderorder=inside['dirlist'].tolist()
	folderindex=[]
	for folder in dirlist:
	  folderindex.append(folderorder.index(folder))
	evaluations.append(evaluation[folderindex])
    if '--electrostatics' in sys.argv:
      if os.path.isfile('Evaluation_'+os.path.splitext(scorename)[0]+'_el_'+ytype+'_'+ycutval+str(probline)+'.npz'):
	inside=np.load('Evaluation_'+os.path.splitext(scorename)[0]+'_el_'+ytype+'_'+ycutval+str(probline)+'.npz')
	evaluation=inside['elevaluation']
	folderorder=inside['dirlist'].tolist()
	folderindex=[]
	for folder in dirlist:
	  folderindex.append(folderorder.index(folder))
	evaluations.append(evaluation[folderindex])
    print '...load done'
  else:
    print 'evaluate ', scorename
    evaluation=[]
    vdwevaluation=[]
    elevaluation=[]
    elon=False
    vdwon=False
    for f,folder in enumerate(dirlist):
      if scorename.split('.')[-1]=='rescore':
	scoring=np.genfromtxt(folder+'/'+scorename)
      elif scorename.split('.')[-1]=='score':
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
	scoring=np.array(scoring)
	if '--vdwpotential' in sys.argv:
	  vdwon=True
	  vdws=np.array(vdws)
	  vsort=np.argsort(vdws)
	  veval=np.cumsum(evals[f][vsort])
	  vdwevaluation.append(veval)
	if '--electrostatics' in sys.argv:
	  elon=True
	  elec=np.array(elec)
	  esort=np.argsort(elec)
	  eleval=np.cumsum(evals[f][esort])
	  elevaluation.append(eleval)
      else:
	print 'datatype ',scorename.split('.')[-1],'not understood'
	sys.exit()
      if len(scoring)==0:
	print scorename,'for ',folder, 'not found'
	sys.exit()
      sort=np.argsort(scoring)
      outeval=np.cumsum(evals[f][sort])
      evaluation.append(outeval)
    evaluation=np.array(evaluation)
    evaluations.append(evaluation)
    np.savez('Evaluation_'+os.path.splitext(scorename)[0]+'_'+ytype+'_'+ycutval+str(probline)+'.npz',dirlist=dirlist,evaluation=evaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    if vdwon:
      vdwevaluation=np.array(vdwevaluation)
      evaluations.append(vdwevaluation)
      np.savez('Evaluation_'+os.path.splitext(scorename)[0]+'_vdw_'+ytype+'_'+ycutval+str(probline)+'.npz',dirlist=dirlist,vdwevaluation=vdwevaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    if elon:
      elevaluation=np.array(elevaluation)
      evaluations.append(elevaluation)
      np.savez('Evaluation_'+os.path.splitext(scorename)[0]+'_el_'+ytype+'_'+ycutval+str(probline)+'.npz',dirlist=dirlist,elevaluation=elevaluation) #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    print '...evaluation done'

if '--deletezeros' in sys.argv:
  eraselist=[]
  for i,folder in enumerate(dirlist):
    erase=ma.nonzero(ma.masked_less(evaluations[0][i],1.))
    if len(erase[0])==0:
      print folder, 'has no good structure'
      eraselist.append(i)

  trainlen=trainlen-len(ma.nonzero(ma.masked_greater_equal(eraselist,trainlen))[0])
  testlen=testlen-len(ma.nonzero(ma.masked_outside(eraselist,len(dirlist)-otherlen-testlen-1,len(dirlist)-otherlen))[0])
  otherlen=otherlen-len(ma.nonzero(ma.masked_less(eraselist,len(dirlist)-otherlen))[0])
  dirlist=np.delete(dirlist,eraselist)
  classes=np.delete(classes, eraselist)
  print np.shape(evaluations)
  evaluations=np.delete(evaluations,eraselist, axis=1)
  print np.shape(evaluations)  
print setname[0], trainlen
print setname[1], testlen
print setname[2], otherlen

if '--plotstructures' in sys.argv:
    x=[]
    simplelist=[]
    class SimpleClass(object):
      pass
    show=[]
    for i in range(len(dirlist)):
      x.append(np.arange(1,len(evaluations[0][i])+1)/float(len(evaluations[0][i])))
      showin=[]
      for c in range(numscores):
	showin.append(evaluations[c][i]/float(max(evaluations[c][i]+[1])))
      show.append(showin)
    x=np.array(x)
    show=np.array(show)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.autoscale(True)
    frame = 0
    ax.set_title(dirlist[frame])
    ax.set_xlabel('n')
    ax.set_ylabel('$n_{pos}/n^{total}_{pos}$')
    ax.set_xscale('log')
    ax.set_yticks(np.arange(0,1.,0.1))
    ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
    ax.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
    ax.grid(which='major', alpha=0.7)
    ax.grid(which='minor', alpha=0.4)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    for count in xrange(numscores):
      lnlist = SimpleClass()
      ln, = ax.step(x[frame], show[frame][count],where='post', color=cm.jet(1.*count/(numscores-1.)),linewidth=1.5, label=legendnames[count][:100])
      lnlist.attr = ln
      simplelist.append(lnlist)
    ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})#,prop={'size':8})
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Frame', 0, len(dirlist)-1, valinit=0,valfmt='%d')
	    # call back function
    def update(val):
	frame = np.floor(sframe.val)
	for count in xrange(numscores):
	  ln=simplelist[count].attr
	  ln.set_xdata(x[frame])
	  ln.set_ydata(show[frame][count])
	ax.set_title(dirlist[frame])
	#ax.autoscale(True)
	#ax.relim()
	#ax.autoscale_view()
	ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
	ax.set_xscale('log')
	ax.set_yticks(np.arange(0,1.,0.1))
	ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
	ax.grid(which='major', alpha=0.7)
	ax.grid(which='minor', alpha=0.4)
	plt.draw()
	    # connect callback to slider   
    sframe.on_changed(update)
    plt.show()
  
if '--atleast' in sys.argv:
    atleast=np.zeros((3,numscores,100000))
    atleastplot=np.zeros((3,numscores,100000))
    atleast1=np.zeros((3,numscores,100000))
    atleast2=np.zeros((3,numscores,100000))    
    if benchon and classon:
      trainlen0 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],0)))
      trainlen1 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],1)))
      trainlen2 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],2)))
      testlen0 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],0)))
      testlen1 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],1)))
      testlen2 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],2)))
    else:
      trainlen0 = trainlen
      testlen0 =testlen
      trainlen1=0
      trainlen2=0
      testlen1=0
      testlen2=0
      
    rank=np.arange(1,100001)
    if '--double%' in sys.argv:
      rank=rank*100./100000.
      for count in range(numscores):
	for c in range(trainlen):
	  sort=np.arange(100000, dtype=int)*len(evaluations[count][c])/100000
	  mask=evaluations[count][c][sort]
	  if len(ma.nonzero(ma.masked_less(mask,1.))[0])!=0:
	    if classes[c]==0:
	      atleast[0,count,ma.nonzero(ma.masked_less(mask,1.))[0][0]]+=1./trainlen0
	    if classes[c]==1:
	      atleast1[0,count,ma.nonzero(ma.masked_less(mask,1.))[0][0]]+=1./trainlen1
	    if classes[c]==2:
	      atleast2[0,count,ma.nonzero(ma.masked_less(mask,1.))[0][0]]+=1./trainlen2
	atleastplot[0,count]=(atleast[0,count]*trainlen0+atleast1[0,count]*trainlen1+atleast2[0,count]*trainlen2)/trainlen
	for c in range(trainlen,trainlen+testlen):
	  sort=np.arange(100000, dtype=int)*len(evaluations[count][c])/100000
	  mask=evaluations[count][c][sort]
	  if len(ma.nonzero(ma.masked_less(mask, 1.))[0])!=0:
	    if classes[c] == 0:
	      atleast[1,count,ma.nonzero(ma.masked_less(mask, 1.))[0][0]]+=1./testlen0
	    if classes[c] == 1:
	      atleast1[1,count,ma.nonzero(ma.masked_less(mask, 1.))[0][0]]+=1./testlen1
	    if classes[c] == 2:
	      atleast2[1,count,ma.nonzero(ma.masked_less(mask, 1.))[0][0]]+=1./testlen2  
	atleastplot[1,count]=(atleast[1,count]*testlen0+atleast1[1,count]*testlen1+atleast2[1,count]*testlen2)/testlen
	for c in range(trainlen+testlen,len(dirlist)):
	  sort=np.arange(100000, dtype=int)*len(evaluations[count][c])/100000
	  mask=evaluations[count][c][sort]
	  if len(ma.nonzero(ma.masked_less(mask, 1.))[0])!=0:
	    atleastplot[2,count,ma.nonzero(ma.masked_less(mask, 1.))[0][0]]+=1./otherlen
    else:
      for count in range(numscores):
	for c in range(trainlen):
	  if len(ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0])!=0:
	    if classes[c] == 0:
	      atleast[0,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./trainlen0
	    if classes[c] == 1:
	      atleast1[0,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./trainlen1
	    if classes[c] == 2:
	      atleast2[0,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./trainlen2
	atleastplot[0,count]=(atleast[0,count]*trainlen0+atleast1[0,count]*trainlen1+atleast2[0,count]*trainlen2)/trainlen
	for c in range(trainlen,trainlen+testlen):
	  if len(ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0])!=0:
	    if classes[c]==0:
	      atleast[1,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./testlen0
	    if classes[c]==1:
	      atleast1[1,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./testlen1
	    if classes[c]==2:
	      atleast2[1,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./testlen2
	atleastplot[1,count]=(atleast[1,count]*testlen0+atleast1[1,count]*testlen1+atleast2[1,count]*testlen2)/testlen
	for c in range(trainlen+testlen, len(dirlist)):
	  if len(ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0])!=0:
	    atleastplot[2,count,ma.nonzero(ma.masked_less(evaluations[count][c],1.))[0][0]]+=1./otherlen		       
    atleast=np.cumsum(atleast,axis=2)
    atleast1=np.cumsum(atleast1,axis=2)
    atleast2=np.cumsum(atleast2,axis=2)
    atleastplot = np.cumsum(atleastplot, axis=2)
    fig=plt.figure(setname[0]+'-atleast_1_in_'+str(trainlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
    ax=fig.add_subplot(111)
    ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
    ax.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$',rotation=0, fontsize = 15)
    ax.yaxis.set_label_coords(-0.12,0.5)
    ax.set_ylim([0,1])
    ax.set_xscale('log') #'symlog', linthreshx=0.1, linscale=1.)
    p=[ax.step(rank,atleastplot[0,i,:nstruc],where='post', color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
    if '--legendinside' in sys.argv:
      ax.legend(loc=0, prop={'size':10})
    else:
      ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    
    if '--double%' in sys.argv:
      ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
      ax.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$',rotation=0, fontsize=15)
      ax.set_yticks(np.arange(0,1.,0.1))
      ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
      ax.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
    ax.grid(which='major', alpha=0.7)
    ax.grid(which='minor', alpha=0.4)
    if benchon == True or classon == True:
      fig2=plt.figure(setname[1]+'-atleast_1_in_'+str(testlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax2=fig2.add_subplot(111)
      ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      ax2.set_ylim([0,1])
      ax2.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
      ax2.yaxis.set_label_coords(-0.12,0.5)
      ax2.set_xscale('log')
      p2=[ax2.step(rank,atleastplot[1,i,:nstruc],where='post', color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	ax2.legend(loc=0, prop={'size':10})
      else:
	ax2.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
      if '--double%' in sys.argv:
	ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	ax2.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
	ax2.yaxis.set_label_coords(-0.12,0.5)
	ax2.set_yticks(np.arange(0,1.,0.1))
	ax2.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax2.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
      ax2.grid(which='major', alpha=0.7)
      ax2.grid(which='minor', alpha=0.4)
    if classon == True and benchon == False:
      fig5=plt.figure(setname[2]+'-atleast_1_in_'+str(otherlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax5=fig5.add_subplot(111)
      ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      ax5.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
      ax5.set_xscale('log')
      p5=[ax5.step(rank,atleastplot[2,i,:nstruc],where='post', color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	ax5.legend(loc=0, prop={'size':10})
      else:
	ax5.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
      if '--double%' in sys.argv:
	ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}/\mathrm{n}_\mathrm{total}$', fontsize=15)
	#ax5.set_ylabel('$n_{complex}^{pos}/n_{complex}^{total}$')
	ax5.set_yticks(np.arange(0,1.,0.1))
	ax5.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax5.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
      ax5.yaxis.set_label_coords(-0.12,0.5)
      ax5.grid(which='major', alpha=0.7)
      ax5.grid(which='minor', alpha=0.4)
    
    
    
    if '--barchart' in sys.argv:
      mean=[]
      meana=[]
      meanb=[]
      print 'percentage of complexes where at least one structures were found'
      fig3=plt.figure(setname[0]+'-atleast_1_in_'+str(trainlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
      ax3=fig3.add_subplot(111)
      ax3.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      ax3.set_ylim([0,1])
      ax3.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
      width=0.9/numscores
      index=np.arange(4)
      if '--double%' in sys.argv:
	print setname[0]
	ax3.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
	ax3.set_xticks(np.arange(0.45,4.45,1.))
	ax3.set_xticklabels( ('top 1','0.1%', '1%', '10%') )
	for count in range(numscores):
	  mean.append([atleast[0,count,0], atleast[0,count,99], atleast[0,count,999],atleast[0,count,9999]])
	  print "%7s" % str(round(100*atleast[0,count,0],2)), "%6s" % str(round(100*atleast[0,count,99],2)), "%6s" % str(round(100*atleast[0,count,999],2)), "%6s" % str(round(100*atleast[0,count,1999],2)), "%6s" % str(round(100*atleast[0,count,4999],2)), "%6s" % str(round(100*atleast[0,count,9999],2)),"%6s" % str(round(100*atleast[0,count,19999],2)), legendnames[count]
	if benchon and classon:
	  print classname[1]
	  for count in range(numscores):
	    meana.append([atleast1[0,count,0], atleast1[0,count,99], atleast1[0,count,999],atleast1[0,count,9999]])
	    print "%7s" % str(round(100*atleast1[0,count,0],2)), "%6s" % str(round(100*atleast1[0,count,99],2)), "%6s" % str(round(100*atleast1[0,count,999],2)), "%6s" % str(round(100*atleast1[0,count,1999],2)), "%6s" % str(round(100*atleast1[0,count,4999],2)), "%6s" % str(round(100*atleast1[0,count,9999],2)),"%6s" % str(round(100*atleast1[0,count,19999],2)), legendnames[count]
	  print classname[2]
	  for count in range(numscores):
	    meanb.append([atleast2[0,count,0], atleast2[0,count,99], atleast2[0,count,999],atleast2[0,count,9999]])
	    print "%7s" % str(round(100*atleast2[0,count,0],2)), "%6s" % str(round(100*atleast2[0,count,99],2)), "%6s" % str(round(100*atleast2[0,count,999],2)), "%6s" % str(round(100*atleast2[0,count,1999],2)), "%6s" % str(round(100*atleast2[0,count,4999],2)), "%6s" % str(round(100*atleast2[0,count,9999],2)),"%6s" % str(round(100*atleast2[0,count,19999],2)), legendnames[count]

      else:
	print setname[0]
	print '| top 1 |  10  |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	ax3.set_xticks(np.arange(0.45,4.45,1.))
	ax3.set_xticklabels( ('1', '10', '100', '1000') )	
	for count in range(numscores):
	  mean.append([atleast[0,count,0], atleast[0,count,9], atleast[0,count,99],atleast[0,count,999]])
	  print "%7s" % str(round(100*atleast[0,count,0],2)), "%6s" % str(round(100*atleast[0,count,9],2)), "%6s" % str(round(100*atleast[0,count,99],2)), "%6s" % str(round(100*atleast[0,count,199],2)), "%6s" % str(round(100*atleast[0,count,499],2)), "%6s" % str(round(100*atleast[0,count,999],2)),"%6s" % str(round(100*atleast[0,count,1999],2)), legendnames[count]
	if benchon and classon:
	  print classname[1]
	  for count in range(numscores):
	    meana.append([atleast1[0,count,0], atleast1[0,count,9], atleast1[0,count,99],atleast1[0,count,999]])
	    print "%7s" % str(round(100*atleast1[0,count,0],2)), "%6s" % str(round(100*atleast1[0,count,9],2)), "%6s" % str(round(100*atleast1[0,count,99],2)), "%6s" % str(round(100*atleast1[0,count,199],2)), "%6s" % str(round(100*atleast1[0,count,499],2)), "%6s" % str(round(100*atleast1[0,count,999],2)),"%6s" % str(round(100*atleast1[0,count,1999],2)), legendnames[count]
	  print classname[2]
	  for count in range(numscores):
	    meanb.append([atleast2[0,count,0], atleast2[0,count,9], atleast2[0,count,99],atleast2[0,count,999]])
	    print "%7s" % str(round(100*atleast2[0,count,0],2)), "%6s" % str(round(100*atleast2[0,count,9],2)), "%6s" % str(round(100*atleast2[0,count,99],2)), "%6s" % str(round(100*atleast2[0,count,199],2)), "%6s" % str(round(100*atleast2[0,count,499],2)), "%6s" % str(round(100*atleast2[0,count,999],2)),"%6s" % str(round(100*atleast2[0,count,1999],2)), legendnames[count]
      if benchon and classon:
	sidewidth =width*0.8/3.
      else:
	sidewidth = width
      ax3.yaxis.set_label_coords(-0.12,0.5)
      p3=[ax3.bar(index+i*width,mean[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      if benchon and classon:
	p3a=[ax3.bar(index+i*width+sidewidth,meana[i], sidewidth, alpha=0.7, color=cm.jet((1.*i)/(numscores-1))) for i in range(numscores)]
	p3b=[ax3.bar(index+i*width+2.*sidewidth,meanb[i], sidewidth, alpha=0.7, color=cm.jet((1.*i)/(numscores-1))) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	ax3.legend(loc=0, prop={'size':10})
      else:
	ax3.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})

      if benchon or classon:    
	mean2=[]
	mean2a=[]
	mean2b=[]
	fig4=plt.figure(setname[1]+'-atleast_1_in_'+str(testlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
	ax4=fig4.add_subplot(111)
	ax4.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
	ax4.set_ylim([0,1])
	ax4.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
	width=0.9/numscores
	index=np.arange(4)
	if '--double%' in sys.argv:
	  print setname[1]
	  ax4.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	  print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
	  ax4.set_xticks(np.arange(0.45,4.45,1.))
	  ax4.set_xticklabels( ('top 1','0.1%', '1%', '10%') )
	  for count in range(numscores):
	    mean2.append([atleast[1,count,0], atleast[1,count,99], atleast[1,count,999], atleast[1,count,9999]])
	    print "%7s" % str(round(100*atleast[1,count,0],2)), "%6s" % str(round(100*atleast[1,count,99],2)), "%6s" % str(round(100*atleast[1,count,999],2)), "%6s" % str(round(100*atleast[1,count,1999],2)), "%6s" % str(round(100*atleast[1,count,4999],2)), "%6s" % str(round(100*atleast[1,count,9999],2)),"%6s" % str(round(100*atleast[1,count,19999],2)), legendnames[count]
	  if benchon and classon:
	    print classname[1]
	    for count in range(numscores):
	      mean2a.append([atleast1[1,count,0], atleast1[1,count,99], atleast1[1,count,999],atleast1[1,count,9999]])
	      print "%7s" % str(round(100*atleast1[1,count,0],2)), "%6s" % str(round(100*atleast1[1,count,99],2)), "%6s" % str(round(100*atleast1[1,count,999],2)), "%6s" % str(round(100*atleast1[1,count,1999],2)), "%6s" % str(round(100*atleast1[1,count,4999],2)), "%6s" % str(round(100*atleast1[1,count,9999],2)),"%6s" % str(round(100*atleast1[1,count,19999],2)), legendnames[count]
	    print classname[2]
	    for count in range(numscores):
	      mean2b.append([atleast2[1,count,0], atleast2[1,count,99], atleast2[1,count,999],atleast2[1,count,9999]])
	      print "%7s" % str(round(100*atleast2[1,count,0],2)), "%6s" % str(round(100*atleast2[1,count,99],2)), "%6s" % str(round(100*atleast2[1,count,999],2)), "%6s" % str(round(100*atleast2[1,count,1999],2)), "%6s" % str(round(100*atleast2[1,count,4999],2)), "%6s" % str(round(100*atleast2[1,count,9999],2)),"%6s" % str(round(100*atleast2[1,count,19999],2)), legendnames[count]

	else:
	  print setname[1]
	  print '| top 1 |  10  |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	  ax4.set_xticks(np.arange(0.45,4.45,1.))
	  ax4.set_xticklabels( ('1', '10', '100', '1000') )	
	  for count in range(numscores):
	    mean2.append([atleast[1,count,0], atleast[1,count,9], atleast[1,count,99],atleast[1,count,999]])
	    print "%7s" % str(round(100*atleast[1,count,0],2)), "%6s" % str(round(100*atleast[1,count,9],2)), "%6s" % str(round(100*atleast[1,count,99],2)), "%6s" % str(round(100*atleast[1,count,199],2)), "%6s" % str(round(100*atleast[1,count,499],2)), "%6s" % str(round(100*atleast[1,count,999],2)),"%6s" % str(round(100*atleast[1,count,1999],2)), legendnames[count]
	  if benchon and classon:
	    print classname[1]
	    for count in range(numscores):
	      mean2a.append([atleast1[1,count,0], atleast1[1,count,9], atleast1[1,count,99],atleast1[1,count,999]])
	      print "%7s" % str(round(100*atleast1[1,count,0],2)), "%6s" % str(round(100*atleast1[1,count,9],2)), "%6s" % str(round(100*atleast1[1,count,99],2)), "%6s" % str(round(100*atleast1[1,count,199],2)), "%6s" % str(round(100*atleast1[1,count,499],2)), "%6s" % str(round(100*atleast1[1,count,999],2)),"%6s" % str(round(100*atleast1[1,count,1999],2)), legendnames[count]
	    print classname[2]
	    for count in range(numscores):
	      mean2b.append([atleast2[1,count,0], atleast2[1,count,9], atleast2[1,count,99],atleast2[1,count,999]])
	      print "%7s" % str(round(100*atleast2[1,count,0],2)), "%6s" % str(round(100*atleast2[1,count,9],2)), "%6s" % str(round(100*atleast2[1,count,99],2)), "%6s" % str(round(100*atleast2[1,count,199],2)), "%6s" % str(round(100*atleast2[1,count,499],2)), "%6s" % str(round(100*atleast2[1,count,999],2)),"%6s" % str(round(100*atleast2[1,count,1999],2)), legendnames[count]
	if benchon and classon:
	  sidewidth =width*0.8/3.
	else:
	  sidewidth = width
	ax4.yaxis.set_label_coords(-0.12,0.5)
	p4=[ax4.bar(index+i*width,mean2[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
	if benchon and classon:
	  p4a=[ax4.bar(index+i*width+sidewidth,mean2a[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
	  p4b=[ax4.bar(index+i*width+2*sidewidth,mean2b[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
	if '--legendinside' in sys.argv:
	  ax4.legend(loc=0, prop={'size':10})
	else:
	  ax4.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
      if classon and benchon == False:    
	mean3=[]
	fig6=plt.figure(setname[2]+'-atleast_1_in_'+str(otherlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
	ax6=fig6.add_subplot(111)
	ax6.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
	ax6.set_ylabel('$\mathrm{n}_\mathrm{c}^\mathrm{pos}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', rotation=0, fontsize=15)
	width=0.9/numscores
	index=np.arange(7)
	if '--double%' in sys.argv:
	  print setname[2]
	  ax6.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	  print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
	  ax6.set_xticklabels( ('top 1', '0.1%', '1%', '2%', '5%', '10%', '20%') )
	  for count in range(numscores):
	    mean3.append([atleast[2,count,0], atleast[2,count,99], atleast[2,count,999],atleast[2,count,1999],atleast[2,count,4999],atleast[2,count,9999],atleast[2,count,19999]])
	    print "%7s" % str(round(100*mean3[count][0],2)), "%6s" % str(round(100*mean3[count][1],2)), "%6s" % str(round(100*mean3[count][2],2)), "%6s" % str(round(100*mean3[count][3],2)), "%6s" % str(round(100*mean3[count][4],2)), "%6s" % str(round(100*mean3[count][5],2)), "%6s" % str(round(100*mean3[count][6],2)), legendnames[count]

	else:
	  print setname[2]
	  print '| top 1 |  10  |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	  ax6.set_xticklabels( ('1', '10', '100', '200', '500', '1000', '2000') )	
	  for count in range(numscores):
	    mean3.append([atleast[2,count,0], atleast[2,count,9], atleast[2,count,99],atleast[2,count,199],atleast[2,count,499],atleast[2,count,999],atleast[2,count,1999]])
	    print "%7s" % str(round(100*mean3[count][0],2)), "%6s" % str(round(100*mean3[count][1],2)), "%6s" % str(round(100*mean3[count][2],2)), "%6s" % str(round(100*mean3[count][3],2)), "%6s" % str(round(100*mean3[count][4],2)), "%6s" % str(round(100*mean3[count][5],2)), "%6s" % str(round(100*mean3[count][6],2)), legendnames[count]
	ax6.yaxis.set_label_coords(-0.12,0.5)
	p6=[ax6.bar(index+i*width,mean3[i], width, alpha=0.6, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
	if '--legendinside' in sys.argv:
	  ax6.legend(loc=0, prop={'size':10})
	else:
	  ax6.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    plt.show()
 
if '--averages' in sys.argv:
    rank=np.arange(1,100001)
    averages=np.zeros((3,numscores,100000))
    averagesplot=np.zeros((3,numscores,100000))
    averages1=np.zeros((3,numscores,100000))
    averages2=np.zeros((3,numscores,100000))    
    if benchon and classon:
      trainlen0 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],0)))
      trainlen1 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],1)))
      trainlen2 = len(ma.compressed(ma.masked_not_equal(classes[:trainlen],2)))
      testlen0 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],0)))
      testlen1 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],1)))
      testlen2 = len(ma.compressed(ma.masked_not_equal(classes[trainlen:],2)))
    else:
      trainlen0 = trainlen
      testlen0 =testlen
      trainlen1=0
      trainlen2=0
      testlen1=0
      testlen2=0
    if '--%' in sys.argv:
      for count in range(numscores):
	for c in range(trainlen):
	  if classes[c]==0:
	    averages[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(trainlen0*max(np.append(evaluations[count][c],[1])))
	  if classes[c]==1:
	    averages1[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(trainlen1*max(np.append(evaluations[count][c],[1])))
	  if classes[c]==2:
	    averages2[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(trainlen2*max(np.append(evaluations[count][c],[1])))
	averagesplot[0,count]=(averages[0,count]*trainlen0+averages1[0,count]*trainlen1+averages2[0,count]*trainlen2)/trainlen    
	for c in range(trainlen,trainlen+testlen):
	  if classes[c]==0:
	    averages[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(testlen0*max(np.append(evaluations[count][c],[1])))
	  if classes[c]==1:
	    averages1[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(testlen1*max(np.append(evaluations[count][c],[1])))
	  if classes[c]==2:
	    averages2[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(testlen2*max(np.append(evaluations[count][c],[1])))
	averagesplot[1,count]=(averages[1,count]*testlen0+averages1[1,count]*testlen1+averages2[1,count]*testlen2)/testlen 	
	for c in range(trainlen+testlen,len(dirlist)):
	  averagesplot[2,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/(otherlen*max(np.append(evaluations[count][c],[1])))
    elif '--double%' in sys.argv:
      rank=rank*100./float(nstruc)
      for count in range(numscores):
	for c in range(trainlen):
	  sort=np.arange(nstruc,dtype=int)*len(evaluations[count][c])/nstruc
	  #sort=sort.astype(int)
	  #print sort
	  if classes[c]==0:
	    averages[0,count]+=evaluations[count][c][sort]/float(trainlen0*max(np.append(evaluations[count][c],[1.])))
	  if classes[c]==1:
	    averages1[0,count]+=evaluations[count][c][sort]/float(trainlen1*max(np.append(evaluations[count][c],[1.])))
	  if classes[c]==2:
	    averages2[0,count]+=evaluations[count][c][sort]/float(trainlen2*max(np.append(evaluations[count][c],[1.])))
	averagesplot[0,count]=(averages[0,count]*trainlen0+averages1[0,count]*trainlen1+averages2[0,count]*trainlen2)/trainlen
	for c in range(trainlen,trainlen+testlen):
	  sort=np.arange(nstruc,dtype= int)*len(evaluations[count][c])/nstruc
	  if classes[c]==0:
	    averages[1,count]+=evaluations[count][c][sort]/float(testlen0*max(np.append(evaluations[count][c],[1.])))
	  if classes[c]==1:
	    averages1[1,count]+=evaluations[count][c][sort]/float(testlen1*max(np.append(evaluations[count][c],[1.])))
	  if classes[c]==2:
	    averages2[1,count]+=evaluations[count][c][sort]/float(testlen2*max(np.append(evaluations[count][c],[1.])))
	averagesplot[1,count]=(averages[1,count]*testlen0+averages1[1,count]*testlen1+averages2[1,count]*testlen2)/testlen 	   
	for c in range(trainlen+testlen, len(dirlist)):
	  sort=np.arange(nstruc,dtype= int)*len(evaluations[count][c])/nstruc
	  averagesplot[2,count]+=evaluations[count][c][sort]/float(otherlen*max(np.append(evaluations[count][c],[1.])))
    else:
      for count in range(numscores):
	for c in range(trainlen):
	  if classes[c]==0:
	    averages[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/trainlen0
	  if classes[c]==1:
	    averages1[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/trainlen1
	  if classes[c]==2:
	    averages2[0,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/trainlen2
	averagesplot[0,count]=(averages[0,count]*trainlen0+averages1[0,count]*trainlen1+averages2[0,count]*trainlen2)/trainlen
	for c in range(trainlen,trainlen+testlen):
	  if classes[c]==0:
	    averages[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/testlen0
	  if classes[c]==1:
	    averages1[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/testlen1
	  if classes[c]==2:
	    averages2[1,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/testlen2
	averagesplot[1,count]=(averages[1,count]*testlen0+averages1[1,count]*testlen1+averages2[1,count]*testlen2)/testlen
	for c in range(trainlen+testlen,len(dirlist)):
	  averagesplot[2,count]+=np.append(evaluations[count][c],max(evaluations[count][c])*np.ones(100000-len(evaluations[count][c])))/otherlen
    
    
    
    fig=plt.figure(setname[0]+'-average_number_goodstructures_in_'+str(trainlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
    ax=fig.add_subplot(111)
    ax.set_xscale('log')
    if '--double%' in sys.argv:
      ax.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
      ax.set_yticks(np.arange(0,1.,0.1))
      ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
      ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]', fontsize = 15)
      ax.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
    elif '--%' in sys.argv:
      ax.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
      ax.set_yticks(np.arange(0,1.,0.1))
      ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
      ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
    else:
      ax.set_ylabel('$E[n_{pos}]$')
      ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
    ax.grid(which='major', alpha=0.7)
    ax.grid(which='minor', alpha=0.4)
    ax.yaxis.set_label_coords(-0.12,0.5)
    p=[ax.step(rank,averages[0,i,:nstruc],where='post', color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
    if '--legendinside' in sys.argv:
      ax.legend(loc=0, prop={'size':10})
    else:
      ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    if benchon or classon:
      fig2=plt.figure(setname[1]+'-average_number_goodstructures_in_'+str(testlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax2=fig2.add_subplot(111)
      ax2.set_xscale('log')
      if '--double%' in sys.argv:
	ax2.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	ax2.set_yticks(np.arange(0,1.,0.1))
	ax2.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	ax2.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
      elif '--%' in sys.argv:
	ax2.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	ax2.set_yticks(np.arange(0,1.,0.1))
	ax2.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      else:
	ax2.set_ylabel('$E[n_{pos}]$')
	ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      ax2.grid(which='major', alpha=0.7)
      ax2.grid(which='minor', alpha=0.4)
      ax2.yaxis.set_label_coords(-0.12,0.5)
      p2=[ax2.step(rank,averages[1,i,:nstruc], where='post',color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	  ax2.legend(loc=0, prop={'size':10})
      else:
	  ax2.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    if classon and benchon == False:
      fig5=plt.figure(setname[2]+'-average_number_goodstructures_in_'+str(otherlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax5=fig5.add_subplot(111)
      ax5.set_xscale('log')
      if '--double%' in sys.argv:
	ax5.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	ax5.set_yticks(np.arange(0,1.,0.1))
	ax5.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	ax5.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
      elif '--%' in sys.argv:
	ax5.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	ax5.set_yticks(np.arange(0,1.,0.1))
	ax5.set_yticks(np.arange(0, 1., 0.05), minor=True)
	ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      else:
	ax5.set_ylabel('E([n_{total}]')
	ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
      ax5.grid(which='major', alpha=0.7)
      ax5.grid(which='minor', alpha=0.4)
      ax5.yaxis.set_label_coords(-0.12,0.5)
      p5=[ax5.step(rank,averages[2,i,:nstruc], where='post',color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	  ax5.legend(loc=0, prop={'size':10})
      else:
	  ax5.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    
    
    
    if '--barchart' in sys.argv:
      mean=[]
      meana=[]
      meanb=[]
      print 'Average percentage of structures found'
      fig3=plt.figure(setname[0]+'-average_number_goodstructures_in_'+str(trainlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
      ax3=fig3.add_subplot(111)
      ax3.set_ylim([0,1])
      width=0.9/numscores
      index=np.arange(4)
      ax3.set_xticks(np.arange(0.45,4.45,1.))
      if '--double%' in sys.argv:
	ax3.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	ax3.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	ax3.set_xticklabels( ('0.1%', '1%', '5%', '10%') )
	print setname[0]
	print '| 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
	for count in range(numscores):
	  mean.append([averages[0,count,99], averages[0,count,999],averages[0,count,4999],averages[0,count,9999]])
	  print "%6s" % str(round(100*averages[0,count,99],2)), "%6s" % str(round(100*averages[0,count,999],2)), "%6s" % str(round(100*averages[0,count,1999],2)), "%6s" % str(round(100*averages[0,count,4999],2)), "%6s" % str(round(100*averages[0,count,9999],2)), "%6s" % str(round(100*averages[0,count,19999],2)), legendnames[count]
	if benchon and classon:
	  print classname[1]
	  for count in range(numscores):
	    meana.append([averages1[0,count,99], averages1[0,count,999],averages1[0,count,4999],averages1[0,count,9999]])
	    print "%6s" % str(round(100*averages1[0,count,99],2)), "%6s" % str(round(100*averages1[0,count,999],2)), "%6s" % str(round(100*averages1[0,count,1999],2)), "%6s" % str(round(100*averages1[0,count,4999],2)), "%6s" % str(round(100*averages1[0,count,9999],2)), "%6s" % str(round(100*averages1[0,count,19999],2)), legendnames[count]
	  print classname[2]
	  for count in range(numscores):
	    meanb.append([averages2[0,count,99], averages2[0,count,999],averages2[0,count,4999],averages2[0,count,9999]])
	    print "%6s" % str(round(100*averages2[0,count,99],2)), "%6s" % str(round(100*averages2[0,count,999],2)), "%6s" % str(round(100*averages2[0,count,1999],2)), "%6s" % str(round(100*averages2[0,count,4999],2)), "%6s" % str(round(100*averages2[0,count,9999],2)), "%6s" % str(round(100*averages2[0,count,19999],2)), legendnames[count]
	  
      else:
	print setname[0]
	print '|  10 |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	ax3.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	ax3.set_ylabel('$E[n_{pos}]$')
	ax3.set_xticklabels( ('10', '100', '500', '1000') )
	for count in range(numscores):
	  mean.append([averages[0,count,9], averages[0,count,99],averages[0,count,499],averages[0,count,999]])
	  print "%6s" % str(round(100*averages[0,count,9],2)), "%6s" % str(round(100*averages[0,count,99],2)), "%6s" % str(round(100*averages[0,count,199],2)), "%6s" % str(round(100*averages[0,count,499],2)), "%6s" % str(round(100*averages[0,count,999],2)), "%6s" % str(round(100*averages[0,count,1999],2)), legendnames[count]
	if benchon and classon:
	  print classname[1]
	  for count in range(numscores):
	    meana.append([averages1[0,count,9], averages1[0,count,99],averages1[0,count,499],averages1[0,count,999]])
	    print "%6s" % str(round(100*averages1[0,count,9],2)), "%6s" % str(round(100*averages1[0,count,99],2)), "%6s" % str(round(100*averages1[0,count,199],2)), "%6s" % str(round(100*averages1[0,count,499],2)), "%6s" % str(round(100*averages1[0,count,999],2)), "%6s" % str(round(100*averages1[0,count,1999],2)), legendnames[count]
	  print classname[2]
	  for count in range(numscores):
	    meanb.append([averages2[0,count,9], averages2[0,count,99],averages2[0,count,499],averages2[0,count,999]])
	    print "%6s" % str(round(100*averages2[0,count,9],2)), "%6s" % str(round(100*averages2[0,count,99],2)), "%6s" % str(round(100*averages2[0,count,199],2)), "%6s" % str(round(100*averages2[0,count,499],2)), "%6s" % str(round(100*averages2[0,count,999],2)), "%6s" % str(round(100*averages2[0,count,1999],2)), legendnames[count]
      
      if benchon and classon:
	sidewith=width/3.
      else:
	sidwidth = width
      ax3.yaxis.set_label_coords(-0.12,0.5)
      p3=[ax3.bar(index+i*width,mean[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      if benchon and classon:
	p3a=[ax3.bar(index+i*width+sidewidth,meana[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
	p3b=[ax3.bar(index+i*width+2*sidewidth,meanb[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
      if '--legendinside' in sys.argv:
	  ax3.legend(loc=0, prop={'size':10})
      else:
	  ax3.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
      
      
      if benchon or classon:    
	mean2=[]
	mean2a=[]
	mean2b=[]
	fig4=plt.figure(setname[1]+'-average_number_goodstructures_in_'+str(testlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
	ax4=fig4.add_subplot(111)
	ax4.set_ylim([0,1])
	width=0.9/numscores
	index=np.arange(4)
	ax4.set_xticks(np.arange(0.45,4.45,1.))
	if '--double%' in sys.argv:	
	  ax4.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	  ax4.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	  ax4.set_xticklabels( ('0.1%', '1%', '5%', '10%') )
	  print setname[1]
	  print '| 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
	  for count in range(numscores):
	    mean2.append([averages[1,count,99], averages[1,count,999],averages[1,count,4999],averages[1,count,9999]])
	    print "%6s" % str(round(100*averages[1,count,99],2)), "%6s" % str(round(100*averages[1,count,999],2)), "%6s" % str(round(100*averages[1,count,1999],2)), "%6s" % str(round(100*averages[1,count,4999],2)), "%6s" % str(round(100*averages[1,count,9999],2)), "%6s" % str(round(100*averages[1,count,19999],2)), legendnames[count]
	  if benchon and classon:
	    print classname[1]
	    for count in range(numscores):
	      mean2a.append([averages1[1,count,99], averages1[1,count,999],averages1[1,count,4999],averages1[1,count,9999]])
	      print "%6s" % str(round(100*averages1[1,count,99],2)), "%6s" % str(round(100*averages1[1,count,999],2)), "%6s" % str(round(100*averages1[1,count,1999],2)), "%6s" % str(round(100*averages1[1,count,4999],2)), "%6s" % str(round(100*averages1[1,count,9999],2)), "%6s" % str(round(100*averages1[1,count,19999],2)), legendnames[count]
	    print classname[2]
	    for count in range(numscores):
	      mean2b.append([averages2[1,count,99], averages2[1,count,999],averages2[1,count,4999],averages2[1,count,9999]])
	      print "%6s" % str(round(100*averages2[1,count,99],2)), "%6s" % str(round(100*averages2[1,count,999],2)), "%6s" % str(round(100*averages2[1,count,1999],2)), "%6s" % str(round(100*averages2[1,count,4999],2)), "%6s" % str(round(100*averages2[1,count,9999],2)), "%6s" % str(round(100*averages2[1,count,19999],2)), legendnames[count]
	else:
	  print setname[1]
	  print '|  10 |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	  ax4.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
	  ax4.set_ylabel('$E[n_{pos}]$')
	  ax4.set_xticklabels( ('10', '100', '500', '1000') )
	  for count in range(numscores):
	    mean2.append([averages[1,count,9], averages[1,count,99],averages[1,count,499],averages[1,count,999]])
	    print "%6s" % str(round(100*averages[1,count,9],2)), "%6s" % str(round(100*averages[1,count,99],2)), "%6s" % str(round(100*averages[1,count,199],2)), "%6s" % str(round(100*averages[1,count,499],2)), "%6s" % str(round(100*averages[1,count,999],2)), "%6s" % str(round(100*averages[1,count,1999],2)), legendnames[count]
	  if benchon and classon:
	    print classname[1]
	    for count in range(numscores):
	      mean2a.append([averages1[1,count,9], averages1[1,count,99],averages1[1,count,499],averages1[1,count,999]])
	      print "%6s" % str(round(100*averages1[1,count,9],2)), "%6s" % str(round(100*averages1[1,count,99],2)), "%6s" % str(round(100*averages1[1,count,199],2)), "%6s" % str(round(100*averages1[1,count,499],2)), "%6s" % str(round(100*averages1[1,count,999],2)), "%6s" % str(round(100*averages1[1,count,1999],2)), legendnames[count]
	    print classname[2]
	    for count in range(numscores):
	      mean2b.append([averages2[1,count,9], averages2[1,count,99],averages2[1,count,499],averages2[1,count,999]])
	      print "%6s" % str(round(100*averages2[1,count,9],2)), "%6s" % str(round(100*averages2[1,count,99],2)), "%6s" % str(round(100*averages2[1,count,199],2)), "%6s" % str(round(100*averages2[1,count,499],2)), "%6s" % str(round(100*averages2[1,count,999],2)), "%6s" % str(round(100*averages2[1,count,1999],2)), legendnames[count]
	ax4.yaxis.set_label_coords(-0.12,0.5)
	p4=[ax4.bar(index+i*width,mean2[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
	if benchon and classon:
	  p4a=[ax4.bar(index+i*width+sidewidth,mean2a[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
	  p4b=[ax4.bar(index+i*width+2*sidewidth,mean2b[i], sidewidth, alpha=0.7, color=cm.jet(1.*i/(numscores-1))) for i in range(numscores)]
	if '--legendinside' in sys.argv:
	  ax4.legend(loc=2, prop={'size':10})
	else:
	  ax4.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
      if classon and benchon == False:    
	mean3=[]
	fig6=plt.figure(setname[2]+'-average_number_goodstructures_in_'+str(otherlen)+'_complexes_'+'_'+ytype+'_'+ycutval+str(probline)+'barchart')
	ax6=fig6.add_subplot(111)
	width=0.9/numscores
	index=np.arange(6)
	if '--double%' in sys.argv:	
	  ax6.set_xlabel('$\mathrm{top}$ $\mathrm{n}$ [%]',fontsize=15)
	  ax6.set_ylabel('$\mathrm{n}_\mathrm{pos}/\mathrm{n}^\mathrm{tot}_\mathrm{pos}$',rotation=0, fontsize=15)
	  ax6.set_xticklabels( ('0.1%', '1%', '2%', '5%', '10%', '20%') )
	  print setname[2]
	  print '| 0.1% |   1% |   2% |   5% |  10% |  20% |  scoring function '
	  for count in range(numscores):
	    mean3.append([averages[2,count,99], averages[2,count,999],averages[2,count,1999],averages[2,count,4999],averages[2,count,9999],averages[2,count,19999]])
	    print "%6s" % str(round(100*mean3[count][0],2)), "%6s" % str(round(100*mean3[count][1],2)), "%6s" % str(round(100*mean3[count][2],2)), "%6s" % str(round(100*mean3[count][3],2)), "%6s" % str(round(100*mean3[count][4],2)), "%6s" % str(round(100*mean3[count][5],2)), legendnames[count]

	else:
	  print setname[2]
	  print '|  10 |  100 |  200 |  500 | 1000 | 2000 | scoring function'
	  ax6.set_xlabel('$\mathrm{top}$ $\mathrm{n}$',fontsize=15)
	  ax6.set_ylabel('$E[n_{pos}]$')
	  ax6.set_xticklabels( ('10', '100', '200', '500', '1000', '2000') )
	  for count in range(numscores):
	    mean3.append([averages[2,count,9], averages[2,count,99],averages[2,count,199],averages[2,count,499],averages[2,count,999],averages[2,count,1999]])
	    print "%6s" % str(round(100*mean3[count][0],2)), "%6s" % str(round(100*mean3[count][1],2)), "%6s" % str(round(100*mean3[count][2],2)), "%6s" % str(round(100*mean3[count][3],2)), "%6s" % str(round(100*mean3[count][4],2)), "%6s" % str(round(100*mean3[count][5],2)), legendnames[count]
	ax6.yaxis.set_label_coords(-0.12,0.5)
	p6=[ax6.bar(index+i*width,mean3[i], width, alpha=0.6, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
	if '--legendinside' in sys.argv:
	  ax6.legend(loc=0, prop={'size':10})
	else:
	  ax6.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    plt.show()
  
if '--ROCstructures' in sys.argv:
    roc1=[]
    roc2=[]
    for count in range(numscores):
      roca=[]
      rocb=[]
      for c in range(len(dirlist)):
	rank=np.arange(1,len(evaluations[count][c])+1,dtype=float)
	badstr=len(evaluations[count][c])-max(evaluations[count][c])
	roca.append(evaluations[count][c]/float(max(np.append(evaluations[count][c],[1.]))))
	rocb.append((rank-evaluations[count][c])/float(badstr))
      roc1.append(roca)
      roc2.append(rocb)
      
    roc=np.array([roc2,roc1])
    simplelist=[]
    class SimpleClass(object):
      pass
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.autoscale(True)
    frame = 0
    ax.set_title(dirlist[frame])
    ax.grid(True)
    ax.set_xlabel('$n_{neg}/n_{neg}^{total}$')
    ax.set_ylabel('$n_{pos}/n_{pos}^{total}$')
    plt.subplots_adjust(left=0.25, bottom=0.25)
    for count in xrange(numscores):
      lnlist = SimpleClass()
      ln, = ax.plot(roc[0,count,frame], roc[1,count,frame], color=cm.jet(1.*count/(numscores-1)), label=legendnames[count][:100])
      lnlist.attr = ln
      simplelist.append(lnlist)
    ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Frame', 0, len(dirlist)-1, valinit=0,valfmt='%d')
	    # call back function
    def update(val):
	frame = np.floor(sframe.val)
	for count in xrange(numscores):
	  ln=simplelist[count].attr
	  ln.set_xdata(roc[0,count,frame])
	  ln.set_ydata(roc[1,count,frame])
	ax.set_title(dirlist[frame])
	ax.grid(True)
	ax.autoscale(True)
	ax.relim()
	ax.autoscale_view()
	ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
	plt.draw()
	    # connect callback to slider   
    sframe.on_changed(update)
    plt.show()
    
if '--ROCcurve' in sys.argv:
    roc=np.zeros((2,3,numscores,100000))
    for count in range(numscores):
      for c in range(trainlen):
	rank=np.arange(1,100001)*100./100000.
	sort=np.arange(100000,dtype=int)*len(evaluations[count][c])/100000
	roc[0,0,count]+=evaluations[count][c][sort]/float(trainlen*max(np.append(evaluations[count][c],[1])))
	roc[1,0,count]+=((np.arange(1,len(evaluations[count][c])+1)-evaluations[count][c])[sort])/float(trainlen*(len(evaluations[count][c])-max(evaluations[count][c])))
      for c in range(trainlen,trainlen+testlen):
	rank=np.arange(1,100001)*100./100000.
	sort=np.arange(100000,dtype=int)*len(evaluations[count][c])/100000
	roc[0,1,count]+=evaluations[count][c][sort]/float(testlen*max(np.append(evaluations[count][c],[1])))
	roc[1,1,count]+=(np.arange(1,len(evaluations[count][c])+1)-evaluations[count][c])[sort]/float(testlen*(len(evaluations[count][c])-max(evaluations[count][c])))
      for c in range(trainlen+testlen,len(dirlist)):
	rank=np.arange(1,100001)*100./100000.
	sort=np.arange(100000,dtype=int)*len(evaluations[count][c])/100000
	roc[0,2,count]+=evaluations[count][c][sort]/float(testlen*max(np.append(evaluations[count][c],[1])))
	roc[1,2,count]+=(np.arange(1,len(evaluations[count][c])+1)-evaluations[count][c])[sort]/float(otherlen*(len(evaluations[count][c])-max(evaluations[count][c])))
    
    fig=plt.figure(setname[0]+'-ROC_'+str(trainlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
    ax=fig.add_subplot(111)    
    ax.set_xlabel('$n_{neg}/n_{neg}^{total}$')
    ax.set_ylabel('$n_{pos}/n_{pos}^{total}$')
    p=[ax.plot(roc[1,0,i], roc[0,0,i], color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
    ax.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    if benchon or classon:
      fig2=plt.figure(setname[1]+'-ROC_'+str(testlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax2=fig2.add_subplot(111)
      ax2.set_xlabel('$n_{neg}/n_{neg}^{total}$')
      ax2.set_ylabel('$n_{pos}/n_{pos}^{total}$')
      p2=[ax2.step(roc[1,1,i],roc[0,1,i], color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      ax2.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    if classon and benchon == False:
      fig3=plt.figure(setname[2]+'-ROC_'+str(otherlen)+'_complexes'+'_'+ytype+'_'+ycutval+str(probline))
      ax3=fig3.add_subplot(111)
      ax3.set_xlabel('$n_{neg}/n_{neg}^{total}$')
      ax3.set_ylabel('$n_{pos}/n_{pos}^{total}$')
      p3=[ax3.step(roc[1,2,i],roc[0,2,i], color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
      ax3.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
    plt.show()
    
