#contacts.py 

import sys, os
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt

#give pregrids with calculated contacts in a given cutoff (bsp gaa-bin7A)
pregrid = sys.argv[1]

dirlist = [x for x in os.listdir('.') if os.path.isdir(x)]

if '--evaluate' in sys.argv:
  Rmsd = sys.argv[sys.argv.index('--evaluate')+1]
  Rmsdend = Rmsd.split('.')[-1]
  rmsdcut = float(sys.argv[sys.argv.index('--evaluate')+2])
else:
  print 'give something to evaluate'
  sys.exit()


goodon=False
if os.path.isfile('Contacts_'+os.path.splitext(pregrid)[0]+'_'+os.path.splitext(Rmsd)[0]+str(rmsdcut)+'.dat'):
  goodon=True

if goodon==False:
  allcounts = []    
  rmsds = []
  print 'read in ...'
  for d, fold in enumerate(dirlist):
    counts = np.load(fold+'/'+pregrid)
    counts = counts.reshape((np.shape(counts)[-2],np.shape(counts)[-1]))
    allcounts.append(counts)
    rmsd = np.genfromtxt(fold+'/'+Rmsd)
    if len(np.shape(rmsd))>1:
      rmsds.append(rmsd[:,-1])
    else:
      rmsds.append(rmsd)
  print '... done'
  
  partypes = np.shape(counts)[-1]
  atypes = int(np.sqrt(0.25+2.*partypes)-0.5)
  
  if '--deletezeros' in sys.argv:
    eraselist=[]
    for i,folder in enumerate(dirlist):
      if Rmsdend[1:5] == 'rmsd':
	erase=ma.nonzero(ma.masked_greater(rmsds[i],rmsdcut))
      else:
	erase=ma.nonzero(ma.masked_less(rmsds[i],rmsdcut))
      if len(erase[0])==0:
	print folder, 'has no good structure'
	eraselist.append(i)
    dirlist=np.delete(dirlist,eraselist)
    rmsds=np.delete(rmsds,eraselist,axis=0)
    allcounts=np.delete(allcounts,eraselist,axis=0)

# goodstuctures
if goodon==True:
  data=np.genfromtxt('Contacts_'+os.path.splitext(pregrid)[0]+'_'+os.path.splitext(Rmsd)[0]+str(rmsdcut)+'.dat',skip_header=1)
  contactmeans=data[:-1,0]
  contactstd=data[:-1,1]
  bcontactmeans=data[:-1,2]
  bcontactstd=data[:-1,3]
  meannumcont=data[-1,0]
  stdnumcont= data[-1,1]
  bmeannumcont=data[-1,2]
  bstdnumcont=data[-1,3]
  print 'mean acceptable structures'
  print meannumcont, stdnumcont
  print 'mean incorrect structures'
  print bmeannumcont, bstdnumcont
  partypes = len(contactmeans)
  atypes = int(np.sqrt(0.25+2.*partypes)-0.5)


else:  
  meancount = np.zeros((len(dirlist),partypes))
  numcontacts = np.zeros(len(dirlist))
  varcontacts =np.zeros(len(dirlist))
  stdcount = np.zeros((len(dirlist),partypes))
  for d in range(len(dirlist)):
      if Rmsdend[1:5] == 'rmsd':
	  maske = ma.nonzero(ma.masked_less_equal(rmsds[d],rmsdcut))
      else:
	  maske = ma.nonzero(ma.masked_greater_equal(rmsds[d]+1.,rmsdcut+1.))

      counts = np.delete(allcounts[d], maske[0], axis = 0)
      numcontacts[d] = np.mean(np.sum(counts,axis=1))
      varcontacts[d] = np.var(np.sum(counts,axis=1))
      counts = counts.astype(np.float32)
      counts = (counts.T/ ma.masked_equal(np.sum(counts,axis = 1),0.).filled(1.)).T
      meancount[d] = np.mean(counts, axis =0)
      stdcount[d]= np.var(counts,axis=0)
      
  contactmeans=np.mean(meancount,axis=0)
  contactstd=np.sqrt(np.mean(stdcount,axis=0))
  
  meannumcont=np.mean(numcontacts,axis=0)
  stdnumcont=np.sqrt(np.mean(varcontacts,axis=0))
  print 'average number of contacts for good structures'

  print contactmeans, contactstd
  print meannumcont, stdnumcont

  # badstuctures
  bmeancount = np.zeros((len(dirlist),partypes))
  bstdcount = np.zeros((len(dirlist),partypes))
  bnumcontacts = np.zeros(len(dirlist))
  bvarcontacts = np.zeros(len(dirlist))
  for d in range(len(dirlist)):
      if Rmsdend[1:5] == 'rmsd':
	  maske = ma.nonzero(ma.masked_greater(rmsds[d],rmsdcut))
      else:
	  maske = ma.nonzero(ma.masked_less(rmsds[d]+1.,rmsdcut+1.))

      counts = np.delete(allcounts[d], maske[0], axis = 0)
      bnumcontacts[d] = np.mean(np.sum(counts,axis=1))
      bvarcontacts[d] = np.var(np.sum(counts,axis=1))
      counts = counts.astype(np.float32)
      counts = (counts.T/ ma.masked_equal(np.sum(counts,axis = 1),0.).filled(1.)).T
      bmeancount[d] = np.mean(counts, axis =0)
      bstdcount[d] = np.var(counts, axis =0)
      
  bcontactmeans=np.mean(bmeancount,axis=0)
  bcontactstd=np.sqrt(np.mean(bstdcount,axis=0))
  bmeannumcont=np.mean(bnumcontacts,axis=0)
  bstdnumcont=np.sqrt(np.mean(bvarcontacts,axis=0))
  print 'average number of contacts for bad structures'
  print bcontactmeans, bcontactstd
  print bmeannumcont, bstdnumcont
  
  
  outarray=np.array([np.append(contactmeans,[meannumcont]), np.append(contactstd,[stdnumcont]),np.append(bcontactmeans,[bmeannumcont]), np.append(bcontactstd,[bstdnumcont])]).T
  np.savetxt('Contacts_'+os.path.splitext(pregrid)[0]+'_'+os.path.splitext(Rmsd)[0]+str(rmsdcut)+'.dat',outarray,header='natmean | natstd | falsemean | falsestd')

if '--nativestructure' in sys.argv:
  natgrid=sys.argv[sys.argv.index('--nativestructure')+1]
  if os.path.isfile('Contacts_'+os.path.splitext(natgrid)[0]+'.dat'):
    datanat=np.genfromtxt('Contacts_'+os.path.splitext(natgrid)[0]+'.dat',skip_header=1)
    natcontactsmean=datanat[:-1,0]
    natcontactsstd=datanat[:-1,1]
    meannumbernat=datanat[-1,0]
    stdnumbernat=datanat[-1,1]
    print 'mean native structures'
    print meannumbernat, stdnumbernat
  else:  
    natcontacts=[]
    for d,folder in enumerate(dirlist):
      counts = np.load(folder+'/'+natgrid)
      counts = counts.reshape((np.shape(counts)[-1]))
      natcontacts.append(counts)

    natcontacts=np.array(natcontacts,dtype=np.float32)
    meannumbernat=np.mean(np.sum(natcontacts,axis=1))
    stdnumbernat=np.std(np.sum(natcontacts,axis=1))
    natcontacts=(natcontacts.T/np.sum(natcontacts,axis=1)).T
    natcontactsmean=np.mean(natcontacts,axis=0)
    natcontactsstd=np.std(natcontacts,axis=0)
    print 'contact distribution of native structures'
    print natcontactsmean, natcontactsstd
    print meannumbernat, stdnumbernat
    
    natoutarray=np.array([np.append(natcontactsmean,[meannumbernat]), np.append(natcontactsmean,[stdnumbernat])]).T
    np.savetxt('Contacts_'+os.path.splitext(natgrid)[0]+'.dat',natoutarray,header='nativemean | nativestd ')


if '--namelist' in sys.argv:
    names = np.genfromtxt(sys.argv[sys.argv.index('--namelist')+1], dtype = str)[:,-1]
else:
    names = np.arange(1,atypes+1,1).astype(str)

if atypes*(atypes-1)/2+atypes == partypes:
  xnames = []
  for n in range(atypes):
    for nn in range(n,atypes):
      xnames.append(names[n]+' / '+names[nn])
else:
  xnames = names

xnames=np.array(xnames, dtype=str)

if '--nativestructure' in sys.argv:
  
  if '--norm' in sys.argv:
    bcontactmeans=ma.masked_equal(bcontactmeans, 0).filled(0.000000001)
    bcontactstd=ma.masked_equal(bcontactstd, 0).filled(0.000000001)
    natcontactsmean=natcontactsmean/bcontactmeans
    contactmeans=contactmeans/bcontactmeans
    natcontactsstd=natcontactsstd/bcontactmeans
    contactstd=contactstd/bcontactmeans
    bcontactstd=bcontactstd/bcontactmeans
    bcontactmeans=bcontactmeans/bcontactmeans
  if '--plotcontact' in sys.argv:
    nativesort=np.argsort(-natcontactsmean)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    width=0.25
    xa=np.arange(0.,partypes,1.)
    ax.set_xticks(xa+0.4)
    ax.set_xticklabels(xnames[nativesort],rotation=90)
    p=ax.bar(xa,natcontactsmean[nativesort], width, color='b')
    p1=ax.bar(xa+width,bcontactmeans[nativesort], width, color='r')
    p2=ax.bar(xa+2*width,contactmeans[nativesort], width, color='g')
    
    plt.show()
    
  if '--plotcontactdist' in sys.argv:
    nativesort=np.argsort(-(contactmeans-bcontactmeans))
    fig=plt.figure('best')
    ax=fig.add_subplot(111)
    width=0.3
    if len(nativesort)<40:
      xa=np.arange(len(nativesort))
      ax.set_xticklabels(xnames[nativesort],rotation=-30, ha='left')
      ax.set_ylabel('$\mathrm{bsA}_\mathrm{i}/<\mathrm{bsA}_\mathrm{i}>$', rotation=0, fontsize=15)
      if '--yerr' in sys.argv:
	p=ax.bar(xa,natcontactsmean[nativesort], width, yerr=natcontactsstd, color='b')
	p2=ax.bar(xa+width,contactmeans[nativesort], width,yerr= contactstd, color='g')
	p3=ax.bar(xa+2*width,bcontactmeans[nativesort], width,yerr= bcontactstd, color='r')
      else:
      	p=ax.bar(xa,natcontactsmean[nativesort], width,color='b')
	p2=ax.bar(xa+width,contactmeans[nativesort], width, color='g')
	p3=ax.bar(xa+2*width,bcontactmeans[nativesort], width, color='r')
    else:
      xa=np.arange(0.,20,1.)
      ax.set_xticklabels(xnames[nativesort][:20],rotation=-30, ha='left')
      ax.set_ylabel('$\mathrm{n}_\mathrm{AB}/<\mathrm{n}_\mathrm{AB}>$', rotation=0, fontsize=15)
      if '--yerr' in sys.argv:
	p=ax.bar(xa,natcontactsmean[nativesort][:20], width,yerr=natcontactsstd[:20], color='b')
	p2=ax.bar(xa+width,contactmeans[nativesort][:20], width, yerr= contactstd[:20], color='g')
	p3=ax.bar(xa+2*width,bcontactmeans[nativesort][:20], width,yerr= bcontactstd[:20], color='r')
      else:
	p=ax.bar(xa,natcontactsmean[nativesort][:20], width, color='b')
	p2=ax.bar(xa+width,contactmeans[nativesort][:20], width, color='g')
	p3=ax.bar(xa+2*width,bcontactmeans[nativesort][:20], width, color='r')

    ax.yaxis.set_label_coords(-0.085,0.5)
    ax.set_xticks(xa)

    if len(nativesort)>40:
      fig1=plt.figure('worst')
      ax1=fig1.add_subplot(111)
      if np.max(np.array([bcontactmeans[nativesort][-20:],contactmeans[nativesort][-20:],natcontactsmean[nativesort][-20:]]))<=1.:
	ax1.set_ylim([0,2])
      width=0.3
      xa=np.arange(0.,20,1.)
      ax1.set_ylabel('$\mathrm{n}_\mathrm{AB}/<\mathrm{n}_\mathrm{AB}>$', rotation=0, fontsize=15)
      ax1.yaxis.set_label_coords(-0.085,0.5)
      ax1.set_xticks(xa+0.3)
      ax1.set_xticklabels(xnames[nativesort][-20:][-np.arange(1,21)],rotation=-30, ha='left')
      if '--yerr' in sys.argv:
	pa=ax1.bar(xa,natcontactsmean[nativesort][-20:][-np.arange(1,21)], width, yerr=natcontactsstd[-20:][-np.arange(1,21)], color='b')
	p2a=ax1.bar(xa+width,contactmeans[nativesort][-20:][-np.arange(1,21)], width, yerr=contactstd[-20:][-np.arange(1,21)], color='g')
	p3a=ax1.bar(xa+2*width,bcontactmeans[nativesort][-20:][-np.arange(1,21)], width, yerr=bcontactstd[-20:][-np.arange(1,21)], color='r')
      else:
	pa=ax1.bar(xa,natcontactsmean[nativesort][-20:][-np.arange(1,21)], width, color='b')
	p2a=ax1.bar(xa+width,contactmeans[nativesort][-20:][-np.arange(1,21)], width, color='g')
	p3a=ax1.bar(xa+2*width,bcontactmeans[nativesort][-20:][-np.arange(1,21)], width, color='r')
    plt.show()

if '--plotparameter' in sys.argv:
  parname = sys.argv[sys.argv.index('--plotparameter')+1]
  parameter = np.genfromtxt(parname,skip_header=1)

  if len(np.shape(parameter))>1:
    atyps = len(parameter[0])
    bins = len(parameter)/atyps
    if bins>1:
      print 'only made for one bin'
      sys.exit()
  
    sortpar = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortnames = []
    ind = 0
    parnorm = np.ones((atyps, atyps))
    if '--normparameter' in sys.argv:
      for n in range(atyps):
	for nn in range(n,atyps):
	  parnorm[n,nn]=bcontactmeans[ind]
	  parnorm[nn,n]=bcontactmeans[ind]
	  ind += 1
      parameter=parameter*parnorm*bmeannumcont
    if '--normparameter2' in sys.argv:
      for n in range(atyps):
	for nn in range(n,atyps):
	  parnorm[n,nn]=abs(contactmeans[ind]-bcontactmeans[ind])
	  parnorm[nn,n]=abs(contactmeans[ind]-bcontactmeans[ind])
	  ind += 1
      parameter=parameter*parnorm*bmeannumcont	
    ind=0
    for n in range(atyps):
      for nn in range(n,atyps):
	sortpar[ind] = parameter[n,nn]
	ind += 1

    sort = np.argsort(sortpar)
    sortpar = sortpar[sort]
    xname = xnames[sort]
    print ' unfavorable contacts               |                  favorable contacts '
    for i in range(30):
      print "%32s" % xname[-i-1], "%3.3f" % sortpar[-i-1],"%32s" % xname[i], "%3.3f" % sortpar[i]
          
    fig0=plt.figure(parname+'-best')
    ax0=fig0.add_subplot(111)
    ax0.set_ylabel('$\mathrm{p}$',rotation=0,fontsize=15)
    ax0.xaxis.set_ticks_position('top')
    ax0.set_xticks(np.arange(0,30))
    ax0.set_xticklabels(xname[:30],rotation=30 , ha='left')
    p0=ax0.bar(np.arange(0,30),sortpar[:30],0.5, color='b', alpha=0.7)
    fig1=plt.figure(parname+'-worst')
    ax1=fig1.add_subplot(111)
    ax1.set_ylabel('$\mathrm{p}$',rotation=0, fontsize=15)
    ax1.set_xticks(np.arange(0,30))
    ax1.set_xticklabels(xname[-30:][-np.arange(1,31)], rotation =-30,ha='left')
    p1=ax1.bar(np.arange(0,30),sortpar[-30:][-np.arange(1,31)],0.5, color='b', alpha=0.7)
 
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    p=ax.imshow(parameter,cmap='spectral',interpolation='none')
    ax.set_yticks(np.arange(0, len(parameter), 1))
    ax.set_xticks(np.arange(0, len(parameter[0]), 1))
    ax.set_xticklabels(names, rotation = 90)
    ax.set_yticklabels(names)
    fig.colorbar(p)
    plt.show()
  
  else:
    if '--normparameter' in sys.argv:
      parameter=parameter*bcontactmeans*bmeannumcont
      sorting=np.argsort(-contactmeans/ma.masked_equal(bcontactmeans,0).filled(0.00000001))
      parameter=parameter[sorting]
      xnames = names[sorting]
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    xb=np.arange(len(parameter))
    p=ax.bar(xb,parameter,0.5, color='b', alpha =0.7)
    ax.set_ylabel('$\mathrm{p}$',rotation=0, fontsize=15)
    ax.set_xticks(np.arange(0, len(parameter), 1))
    ax.set_xticklabels(xnames, rotation = -30, ha='left')

    plt.show()
  #if '--native-false' in sys.argv:
    
    
  #if '--near-false' in sys.argv:


if '--plotparametercontribution' in sys.argv:
  parname = sys.argv[sys.argv.index('--plotparametercontribution')+1]
  parameter = np.genfromtxt(parname,skip_header=1)

  if len(np.shape(parameter))>1:
    atyps = len(parameter[0])
    bins = len(parameter)/atyps
    if bins>1:
      print 'only made for one bin'
      sys.exit()
  
    sortpar = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortpar0 = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortpar1 = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortpar2 = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortpar3 = np.zeros(((atyps*(atyps-1)/2)+atyps))
    sortnames = []
    ind = 0
    parnorm = np.ones((atyps, atyps))
    parnorm0 = np.ones((atyps, atyps))
    parnorm1 = np.ones((atyps, atyps))
    parnorm2 = np.ones((atyps, atyps))
    parnorm3 = np.ones((atyps, atyps))
      
    if '--normparameter' in sys.argv:
      for n in range(atyps):
	for nn in range(n,atyps):
	  parnorm1[n,nn]=abs(contactmeans[ind]-bcontactmeans[ind])
	  parnorm1[nn,n]=abs(contactmeans[ind]-bcontactmeans[ind])
	  parnorm0[n,nn]=bcontactmeans[ind]
	  parnorm0[nn,n]=bcontactmeans[ind]
	  parnorm2[n,nn]=contactmeans[ind]
	  parnorm2[nn,n]=contactmeans[ind]
	  parnorm3[n,nn]=natcontactsmean[ind]
	  parnorm3[nn,n]=natcontactsmean[ind]
	  ind += 1
      parameter1=parameter*parnorm1*bmeannumcont
      parameter0=parameter*parnorm0*bmeannumcont
      parameter2=parameter*parnorm2*bmeannumcont
      parameter3=parameter*parnorm3*bmeannumcont
    
    ind=0
    for n in range(atyps):
      for nn in range(n,atyps):
	sortpar0[ind] = parameter0[n,nn]
	sortpar1[ind] = parameter1[n,nn]
	sortpar2[ind] = parameter2[n,nn]
	sortpar3[ind] = parameter3[n,nn]
	ind += 1

    sort = np.argsort(sortpar1)
    sortpar0 = sortpar0[sort]
    sortpar2 = sortpar2[sort]
    sortpar3 = sortpar3[sort]
    xname = xnames[sort]
    print ' unfavorable contacts               |                  favorable contacts '
    for i in range(30):
      print "%32s" % xname[-i-1], "%3.3f" % sortpar[-i-1],"%32s" % xname[i], "%3.3f" % sortpar[i]
          
    fig0=plt.figure(parname+'-best')
    ax0=fig0.add_subplot(111)
    ax0.set_ylabel('$\mathrm{p}$',rotation=0,fontsize=15)
    ax0.xaxis.set_ticks_position('top')
    ax0.set_xticks(np.arange(0,20))
    ax0.set_xticklabels(xname[:20],rotation=30 , ha='left')
    p0=ax0.bar(np.arange(0,20)+0.6,sortpar0[:20],0.3, color='r', alpha=0.7)
    p0b=ax0.bar(np.arange(0,20)+0.3,sortpar2[:20],0.3, color='g', alpha=0.7)
    p0c=ax0.bar(np.arange(0,20),sortpar3[:20],0.3, color='b', alpha=0.7)
    fig1=plt.figure(parname+'-worst')
    ax1=fig1.add_subplot(111)
    ax1.set_ylabel('$\mathrm{p}$',rotation=0, fontsize=15)
    ax1.set_xticks(np.arange(0,20))
    ax1.set_xticklabels(xname[-20:][-np.arange(1,21)], rotation =-30,ha='left')
    p1=ax1.bar(np.arange(0,20)+0.6,sortpar0[-20:][-np.arange(1,21)],0.3, color='r', alpha=0.7)
    p1b=ax1.bar(np.arange(0,20)+0.3,sortpar2[-20:][-np.arange(1,21)],0.3, color='g', alpha=0.7)
    p1b=ax1.bar(np.arange(0,20),sortpar3[-20:][-np.arange(1,21)],0.3, color='b', alpha=0.7)
    
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    p=ax.imshow(parameter1,cmap='spectral',interpolation='none')
    ax.set_yticks(np.arange(0, len(parameter), 1))
    ax.set_xticks(np.arange(0, len(parameter[0]), 1))
    ax.set_xticklabels(names, rotation = 90)
    ax.set_yticklabels(names)
    fig.colorbar(p)
    plt.show() 
    
  
  else:
    if '--normparameter' in sys.argv:
      parameter3=parameter*abs(contactmeans-bcontactmeans)*bmeannumcont
      parameter0=parameter*contactmeans*bmeannumcont
      parameter1=parameter*bcontactmeans*bmeannumcont
      parameter2=parameter*natcontactsmean*bmeannumcont
      sorting=ma.compressed(ma.array(np.argsort(parameter3),mask=ma.getmask(ma.masked_equal(np.sort(parameter3),0))))      
      parameter0=parameter0[sorting]
      parameter1=parameter1[sorting]
      parameter2=parameter2[sorting]
      xnames = names[sorting]
    fig=plt.figure(parname)
    ax=fig.add_subplot(111)
    xb=np.arange(len(parameter0))
    p=ax.bar(xb+0.3,parameter0,0.3, color='g', alpha =0.7)
    p1=ax.bar(xb+0.6,parameter1,0.3, color='r', alpha =0.7)
    p2=ax.bar(xb,parameter2,0.3, color='b', alpha =0.7)
    
    ax.set_ylabel('$\mathrm{p}$',rotation=0, fontsize=15)
    ax.set_xticks(np.arange(0, len(parameter), 1))
    ax.set_xticklabels(xnames, rotation = -30, ha='left')

    plt.show()    
    

