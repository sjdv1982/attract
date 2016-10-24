#compare-nativescores.py

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


nstruc = 100000
evaluations = []
numstrucs = []
for n, scorename in enumerate(scorenames):
  elon=True
  vdwon=True
  if scorename.split('.')[-1]=='dat':
    if '--electrostatics' in sys.argv:
      elon=False
      if os.path.isfile('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_el.dat'):
	elon=True
    if '--vdwpotential' in sys.argv:
      vdwon=False
      if os.path.isfile('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_vdw.dat'):
	vdwon=True	
  if vdwon and elon and os.path.isfile('Native-Evaluation_'+os.path.splitext(scorename)[0]+'.dat'):
    print 'load Native-Evaluation_'+os.path.splitext(scorename)[0]+'.dat ...'
    inside=np.genfromtxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'.dat', skip_header =1, dtype = '|5S')
    folderorder=list(inside[:,0])
    evaluation=inside[:,1].astype(np.int)
    numstruc = inside[:,2].astype(np.int)
    folderindex=[]
    for folder in dirlist:
      folderindex.append(folderorder.index(folder))
    evaluations.append(evaluation[folderindex])
    numstrucs.append(numstruc[folderindex])
    if '--vdwpotential' in sys.argv:
      if os.path.isfile('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_vdw.dat'):
	inside=np.genfromtxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_vdw.dat', skip_header = 1, dtype = '|5S')
	folderorder=list(inside[:,0])
	evaluation=inside[:,1].astype(np.int)
	numstruc = inside[:,2].astype(np.int)
	folderindex=[]
	for folder in dirlist:
	  folderindex.append(folderorder.index(folder))
	evaluations.append(evaluation[folderindex])
	numstrucs.append(numstruc[folderindex])
    if '--electrostatics' in sys.argv:
      if os.path.isfile('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_el.dat'):
	inside=np.genfromtxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_el.dat', skip_header = 1, dtype = '|5S')
	folderorder=list(inside[:,0])
	evaluation=inside[:,1].astype(np.int)
	numstruc = inside[:,2].astype(np.int)
	folderindex=[]
	for folder in dirlist:
	  folderindex.append(folderorder.index(folder))
	evaluations.append(evaluation[folderindex])
	numstrucs.append(numstruc[folderindex])
    print '...load done'
  else:
    print 'evaluate ', scorename
    outarray = np.ndarray((len(dirlist),3),dtype = object)
    outarray[:,0] = dirlist
    evaluation=[]
    vdwevaluation=[]
    elevaluation=[]
    numstruc = []
    elon=False
    vdwon=False
    for f,folder in enumerate(dirlist):
      if scorename.split('.')[-1]=='rescore':
	scoring=np.genfromtxt(folder+'/'+scorename)
	natscore = np.genfromtxt(folder+'/'+os.path.splitext(scorename)[0]+'-native.rescore')
	scoring = np.append(scoring,natscore)
      elif scorename.split('.')[-1]=='score':
	scoring=np.genfromtxt(folder+'/'+scorename)
	natscore = np.genfromtxt(folder+'/'+os.path.splitext(scorename)[0]+'-native.rescore')
	scoring = np.append(scoring,natscore)
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
	natobj = open(folder+'/'+os.path.splitext(scorename)[0]+'-native.rescore', 'r')
	nlines = natobj.readlines()
	for n, nline in enumerate(nlines):
	  if nline[:8] == ' Energy:':
	    nline=nline.strip()
	    nline=nline.split()
	    scoring.append(float(nline[1]))
	    if '--electrostatics' in sys.argv:
	      ele=nlines[n+1].strip()
	      ele=ele.split()
	      elec=np.append(elec,[float(ele[1])])
	    if '--vdwpotential' in sys.argv:
	      v=nlines[n+1].strip()
	      v=v.split()
	      if v[0][0]=='*':
		vdws=np.append(vdws,[10000000])
	      else:
		vdws=np.append(vdws,[float(v[0])])
	scoring=np.array(scoring)
	if '--vdwpotential' in sys.argv:
	  vdwon=True
	  vdws=np.array(vdws)
	  vsort=np.searchsorted(vdws[:-1],vdws[-1])
	  vdwevaluation.append(vsort)
	if '--electrostatics' in sys.argv:
	  elon=True
	  elec=np.array(elec)
	  esort=np.searchsorted(elec[:-1],elec[-1])
	  elevaluation.append(esort)
      else:
	print 'datatype ',scorename.split('.')[-1],'not understood'
	sys.exit()
      if len(scoring)==0:
	print scorename,'for ',folder, 'not found'
	sys.exit()
      #insert native scores
      
      sort=np.searchsorted(scoring[:-1], scoring[-1])
      evaluation.append(sort)
      numstruc.append(len(scoring))
    evaluation=np.array(evaluation)
    evaluations.append(evaluation)
    numstrucs.append(numstruc)
    outarray[:,1] = evaluation
    outarray[:,2] = numstruc
    np.savetxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'.dat',outarray, fmt = ["%5s", "%d", "%d"], header = 'Complex | natrank | numstruc') #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    if vdwon:
      vdwevaluation=np.array(vdwevaluation)
      outarray[:,1] = vdwevaluation
      evaluations.append(vdwevaluation)
      np.savetxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_vdw.dat',outarray, fmt = ["%5s", "%d", "%d"], header = 'Complex | natrank | numstruc') #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    if elon:
      elevaluation=np.array(elevaluation)
      evaluations.append(elevaluation)
      outarray[:,1] = elevaluation
      np.savetxt('Native-Evaluation_'+os.path.splitext(scorename)[0]+'_el.dat',outarray, fmt = ["%5s", "%d", "%d"], header = 'Complex | natrank | numstruc') #np.append(np.reshape(setlength,(len(setlength),1)),evaluation,axis=1))#,comments='', header=head)
    print '...evaluation done'
    
evaluations = np.array(evaluations, dtype = np.float)
numstrucs = np.array(numstrucs, dtype = np.float)

natives=np.zeros((3,numscores,100000))
rank = np.arange(1,100001,1)
if '--double%' in sys.argv:
  for count in range(numscores):
    for c in range(trainlen):
      sort=int(evaluations[count][ c]*100000./numstrucs[count][c])
      natives[0,count,sort:] += 1./trainlen
    for c in range(trainlen,trainlen+testlen):
      sort=int(evaluations[count][ c]*100000./numstrucs[count][c])
      natives[1,count,sort:] += 1./testlen
    for c in range(trainlen+testlen,len(dirlist)):
      sort=int(evaluations[count][ c]*100000./numstrucs[count][c])
      natives[2,count,sort:] += 1./otherlen
else:
  for count in range(numscores):
    for c in range(trainlen):
      sort=int(evaluations[count][ c])
      natives[0,count,sort:] += 1./trainlen
    for c in range(trainlen,trainlen+testlen):
      sort=int(evaluations[count][ c])
      natives[1,count,sort:] += 1./testlen
    for c in range(trainlen+testlen,len(dirlist)):
      sort=int(evaluations[count][ c])
      natives[2,count,sort:] += 1./otherlen	       
    
    
fig=plt.figure(setname[0]+'_nativescore_'+str(trainlen)+'_complexes')
ax=fig.add_subplot(111)
ax.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
ax.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
ax.yaxis.set_label_coords(-0.12,0.5)
ax.set_xscale('log') #'symlog', linthreshx=0.1, linscale=1.)
p=[ax.step(rank,natives[0,i,:nstruc],where='post',  color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
ax.legend(loc=4, prop={'size':10})
ax.set_ylim([0.,1.1])
if '--double%' in sys.argv:
  ax.set_yticks(np.arange(0,1.,0.1))
  ax.set_yticks(np.arange(0, 1., 0.05), minor=True)
  ax.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
ax.grid(which='major', alpha=0.7)
ax.grid(which='minor', alpha=0.4)
if benchon == True or classon == True:
  fig2=plt.figure(setname[1]+'_nativescore_'+str(testlen)+'_complexes')
  ax2=fig2.add_subplot(111)
  ax2.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
  ax2.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
  ax2.yaxis.set_label_coords(-0.12,0.5)
  ax2.set_xscale('log')
  p2=[ax2.step(rank,natives[1,i,:nstruc],where='post',  color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
  ax2.legend(loc=4, prop={'size':10})
  ax2.set_ylim([0.,1.1])
  if '--double%' in sys.argv:
    ax2.set_yticks(np.arange(0,1.,0.1))
    ax2.set_yticks(np.arange(0, 1., 0.05), minor=True)
    ax2.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
  ax2.grid(which='major', alpha=0.7)
  ax2.grid(which='minor', alpha=0.4)
if classon == True:
  fig5=plt.figure(setname[2]+'_nativescore_'+str(testlen)+'_complexes')
  ax5=fig5.add_subplot(111)
  ax5.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
  ax5.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
  ax5.yaxis.set_label_coords(-0.12,0.5)
  ax5.set_xscale('log')
  p5=[ax5.step(rank,natives[2,i,:nstruc],where='post',  color=cm.jet(1.*i/(numscores-1.)), label=legendnames[i]) for i in range(numscores)]
  ax5.legend(loc=4, prop={'size':10})
  ax5.set_ylim([0.,1.1])
  if '--double%' in sys.argv:
    ax5.set_yticks(np.arange(0,1.,0.1))
    ax5.set_yticks(np.arange(0, 1., 0.05), minor=True)
    ax5.set_xticklabels(('0.000001','0.00001','0.0001','0.001', '0.01', '0.1', '1'))
  ax5.grid(which='major', alpha=0.7)
  ax5.grid(which='minor', alpha=0.4)



if '--barchart' in sys.argv:
  mean=[]
  print 'percentage of complexes where native structures were found'
  fig3=plt.figure(setname[0]+'_nativescore_'+str(trainlen)+'_complexes'+'barchart')
  ax3=fig3.add_subplot(111)
  ax3.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
  ax3.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
  ax3.yaxis.set_label_coords(-0.12,0.5)
  ax3.set_ylim([0.,1.1])
  width=0.9/numscores
  index=np.arange(4)
  ax3.set_xticks(np.arange(0.45,4.45,1.))
  if '--double%' in sys.argv:
    print setname[0]
    print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
    ax3.set_xticklabels( ('top 1', '0.1%', '1%', '10%') )
    for count in range(numscores):
      mean.append([natives[0,count,0], natives[0,count,99], natives[0,count,999],natives[0,count,9999]])
      print "%6s" % str(round(100*natives[0,count,0],2)), "%6s" % str(round(100*natives[0,count,99],2)), "%6s" % str(round(100*natives[0,count,999],2)), "%6s" % str(round(100*natives[0,count,1999],2)), "%6s" % str(round(100*natives[0,count,4999],2)), "%6s" % str(round(100*natives[0,count,9999],2)), "%6s" % str(round(100*natives[0,count,19999],2)), legendnames[count]
  else:
    print setname[0]
    print '| top 1 |  10 |  100 |  200 | 500 | 1000 | 2000 | scoring function'
    ax3.set_xticklabels( ('1', '10', '100', '1000') )	
    for count in range(numscores):
      print "%6s" % str(round(100*natives[0,count,0],2)), "%6s" % str(round(100*natives[0,count,9],2)), "%6s" % str(round(100*natives[0,count,99],2)), "%6s" % str(round(100*natives[0,count,199],2)), "%6s" % str(round(100*natives[0,count,499],2)), "%6s" % str(round(100*natives[0,count,999],2)), "%6s" % str(round(100*natives[0,count,1999],2)), legendnames[count]
      mean.append([natives[0,count,0],natives[0,count,9], natives[0,count,99],natives[0,count,999]])
  p3=[ax3.bar(index+i*width,mean[i], width, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
  ax3.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
  if benchon or classon:    
    mean2=[]
    fig4=plt.figure(setname[1]+'_nativescore_'+str(testlen)+'_complexes'+'barchart')
    ax4=fig4.add_subplot(111)
    ax4.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
    ax4.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
    ax4.yaxis.set_label_coords(-0.12,0.5)
    ax4.set_ylim([0.,1.1])
    width=0.9/numscores
    index=np.arange(4)
    ax4.set_xticks(np.arange(0.45,4.45,1.))
    if '--double%' in sys.argv:
      print setname[1]
      print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
      ax4.set_xticklabels( ('top 1', '0.1%', '1%', '10%') )
      for count in range(numscores):
	mean2.append([natives[1,count,0], natives[1,count,99], natives[1,count,999],natives[1,count,9999]])
	print "%6s" % str(round(100*natives[1,count,0],2)), "%6s" % str(round(100*natives[1,count,99],2)), "%6s" % str(round(100*natives[1,count,999],2)), "%6s" % str(round(100*natives[1,count,1999],2)), "%6s" % str(round(100*natives[1,count,4999],2)), "%6s" % str(round(100*natives[1,count,9999],2)), "%6s" % str(round(100*natives[1,count,19999],2)), legendnames[count]
    else:
      print setname[1]
      print '| top 1 |  10 |  100 |  200 | 500 | 1000 | 2000 | scoring function'
      ax4.set_xticklabels( ('1', '10', '100', '1000') )	
      for count in range(numscores):
	print "%6s" % str(round(100*natives[1,count,0],2)), "%6s" % str(round(100*natives[1,count,9],2)), "%6s" % str(round(100*natives[1,count,99],2)), "%6s" % str(round(100*natives[1,count,199],2)), "%6s" % str(round(100*natives[1,count,499],2)), "%6s" % str(round(100*natives[1,count,999],2)), "%6s" % str(round(100*natives[1,count,1999],2)), legendnames[count]
	mean2.append([natives[1,count,0],natives[1,count,9], natives[1,count,99],natives[1,count,999]])
    p4=[ax4.bar(index+i*width,mean2[i], width, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
    ax4.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
  if classon:    
    mean3=[]
    fig6=plt.figure(setname[2]+'_nativescore_'+str(otherlen)+'_complexes'+'barchart')
    ax6=fig6.add_subplot(111)
    ax6.set_xlabel('$\mathrm{top}$ $\mathrm{n}_\mathrm{decoy}$', fontsize=15)
    ax6.set_ylabel('$\mathrm{n}_\mathrm{c}/\mathrm{n}_\mathrm{c}^\mathrm{tot}$', fontsize=15, rotation=0)
    ax6.yaxis.set_label_coords(-0.12,0.5)
    ax6.set_ylim([0.,1.1])
    width=0.9/numscores
    index=np.arange(4)
    ax6.set_xticks(np.arange(0.45,4.45,1.))
    if '--double%' in sys.argv:
      print setname[2]
      print '| top 1 | 0.1% |   1% |   2% |   5% |  10% |  20% | scoring function'
      ax6.set_xticklabels( ('top 1', '0.1%', '1%', '10%') )
      for count in range(numscores):
	mean3.append([natives[2,count,0], natives[2,count,99], natives[2,count,999],natives[2,count,9999]])
	print "%6s" % str(round(100*natives[2,count,0],2)), "%6s" % str(round(100*natives[2,count,99],2)), "%6s" % str(round(100*natives[2,count,999],2)), "%6s" % str(round(100*natives[2,count,1999],2)), "%6s" % str(round(100*natives[2,count,4999],2)), "%6s" % str(round(100*natives[2,count,9999],2)), "%6s" % str(round(100*natives[2,count,19999],2)), legendnames[count]
    else:
      print setname[2]
      print '| top 1 |  10 |  100 |  200 | 500 | 1000 | 2000 | scoring function'
      ax6.set_xticklabels( ('1', '10', '100', '1000') )	
      for count in range(numscores):
	print "%6s" % str(round(100*natives[2,count,0],2)), "%6s" % str(round(100*natives[2,count,9],2)), "%6s" % str(round(100*natives[2,count,99],2)), "%6s" % str(round(100*natives[2,count,199],2)), "%6s" % str(round(100*natives[2,count,499],2)), "%6s" % str(round(100*natives[2,count,999],2)), "%6s" % str(round(100*natives[2,count,1999],2)), legendnames[count]
	mean3.append([natives[2,count,0],natives[2,count,9], natives[2,count,99],natives[2,count,999]])
    p6=[ax6.bar(index+i*width,mean3[i], width, alpha=0.7, color=cm.jet(1.*i/(numscores-1)), label=legendnames[i]) for i in range(numscores)]
    ax6.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
plt.show()
    
    
