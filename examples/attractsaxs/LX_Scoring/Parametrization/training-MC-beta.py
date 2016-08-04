import sys,os
import numpy as np
import numpy.ma as ma
import random as random
import matplotlib.pyplot as plt
import duplicatelib
import fgenlib
np.set_printoptions(precision=3,suppress=True,linewidth=2000)

#run it by: 

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
  elif gridtype == 'interpolate':
    start = float(sys.argv[sys.argv.index('--grid')+4])
    rcutoff = float(sys.argv[sys.argv.index('--grid')+5])
    stepsize = float(sys.argv[sys.argv.index('--grid')+6])
    potsteps = int((rcutoff-start)/stepsize)+1
    Header += ' start: '+str(start)+' rcut: '+str(rcutoff)+' stepsize: '+str(stepsize)
  
else:
  print 'please give a grid to train on'
  sys.exit()

if '--mcparams' in sys.argv:
  stepping = sys.argv[sys.argv.index('--mcparams')+1]
  parstepsize = float(sys.argv[sys.argv.index('--mcparams')+2])
  ansteps = sys.argv[sys.argv.index('--mcparams')+3]
  target=sys.argv[sys.argv.index('--mcparams')+4]
  Header += ' steps: '+stepping+' '+str(parstepsize)+' '+ansteps+' Targetfunction: '+target
  tind = sys.argv.index('--mcparams')+4
  if target == 'rank' or target[:6] == 'refine' or target == 'positiontop':
    bestof = int(sys.argv[tind+1])
    tind +=1
    Header += str(bestof)
  annealing=sys.argv[tind + 1]
  temp=sys.argv[tind + 2]
  mctype=sys.argv[tind + 3]
  Header += ' Annealing: '+annealing+' T0: '+temp+' mctpye: '+mctype
  if mctype == 'interpolate':    
    fgentype = int(sys.argv[tind + 4])
    power1 = float(sys.argv[tind + 5])
    power2 = float(sys.argv[tind + 6])
    Header += ' fgentype: '+str(fgentype)+' power: '+str(power1)+' '+str(power2)
  elif mctype == 'normal' or mctype == 'keepsign':
    changing = sys.argv[tind+4]
    Header += ' changing: '+changing
    if changing == 'saddle':
      power1 = float(sys.argv[tind+5])
      power2 = float(sys.argv[tind+6])
      Header +=' power: '+str(power1)+' '+str(power2)
else:
  print 'insert parameter for montecarlo annealing' 
  sys.exit()
  
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
  
if '--cutoffweight' in sys.argv:
  evcut = float(sys.argv[sys.argv.index('--cutoffweight')+1])

if '--maxweight' in sys.argv:
  maxcut = float(sys.argv[sys.argv.index('--maxweight')+1])
  
if '--preparameter' in sys.argv:
  parmw=sys.argv[sys.argv.index('--preparameter')+1]  # random, sattlepoint, ABC
  bins = int(sys.argv[sys.argv.index('--preparameter')+2])
  paratyps = int(sys.argv[sys.argv.index('--preparameter')+3])
  Header += ' Parameter: '+parmw+' bins: '+str(bins)+' #parameter: '+str(paratyps)
  if parmw == 'saddlepoint' or parmw == 'ABC':
    parname = sys.argv[sys.argv.index('--preparameter')+4]
    Header += parname
else:
  print 'insert paramter file to start or random.bin.atomtypes'
  sys.exit()

constpar = False
if '--searchrange' in sys.argv:
  constpar = True
  parranges = np.ones((bins,2))
  Header += 'constrains: '
  for i in range(bins):
    parranges[i,0]=float(sys.argv[sys.argv.index('--searchrange')+1+i*2])
    parranges[i,1]=float(sys.argv[sys.argv.index('--searchrange')+2+i*2])
    Header += str(parranges[i,0])+' '+str(parranges[i,1])
    
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
  print 'insert name and form for outputfiles'
  sys.exit()
  
cmplx=np.genfromtxt(Complexes, dtype=str)
numcpx=len(cmplx)
structures=int(strux)		#number of decoys taken into account

rmsd=np.zeros((numcpx,structures),dtype=np.float32)		#matrix of all rmsds of the complexes
for d,name in enumerate(cmplx):					#for every complex one makes a matrix with all the counts of atompairs
    #read rmsds
    rmsds = np.genfromtxt(name+'/'+Rmsd)
    if len(np.shape(rmsds))>1:
	rmsd[d]= rmsds[:structures,-1]
    elif len(np.shape(rmsds))==1:
	rmsd[d]= rmsds[:structures]
    if natives:
      if Rmsd.split('.')[-1][1:5]=='rmsd':
	rmsd[d,-1]=0.001
      elif Rmsd.split('.')[-1][:8]=='capstars':
	rmsd[d,-1]==4
      elif Rmsd.split('.')[-1][:4]=='fnat':
	rmsd[d,-1]==1.
  

allcounts = []
weight=np.zeros((numcpx,structures))

for d, name in enumerate(cmplx):
  if gridtype == 'interpolate':
    counts = np.load(name+'/'+Counts, mmap_mode='r')
    counts = counts[:,:structures,:]
  elif gridtype == 'distances':
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
      counts = counts.reshape((bins,structures,paratyps))
    elif len(np.shape(counts)) == 3:
      counts = counts[:,:structures,:]
    else:
      print 'shape of grids not understood', np.shape(counts)
      sys.exit()
  
  if natives:
    if gridtype == 'interpolate':
      natcounts = np.load(name+'/'+nativegrid, mmap_mode='r')
    elif gridtype == 'distances':
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
	natcounts = natcounts.reshape((bins,structures,paratyps))
    counts=np.append(counts[:,:-1,:],natcounts, axis=1)
    
  if dup == 'duplication':
    counts,weight[d]=duplicatelib.duplicate(rmsd[d,:], counts, valuate, prorow)
  
  allcounts.append(counts)

allcounts=np.array(allcounts, dtype = np.float32)

if '--eraseatomtype' in sys.argv:
    atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)
    hydindex=0
    hydindexes=[]
    for n in range(atyps):
	for nn in range(n,atyps):
	    if nn==eratomtype-1:
		hydindexes.append(hydindex)
	    hydindex+=1
    allcounts=np.delete(allcounts, hydindexes, -1)
    paratyps=np.shape(allcounts)[-1]
    atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)

if parmw == 'random':
  potpars = np.random.rand(bins,paratyps)
  if gridtype == 'interpolate':
    potpars[-1] = np.random.choice([-1,1],paratyps)
    potpars[1] = potpars[1]*20.*np.random.random()
    if mctype == 'interpolate':
      pot = fgenlib.fgen(stepsize, start, rcutoff, potsteps, fgentype, potpars, power1, power2, paratyps, -1)
    else: 
      print 'gridtype interpolate works only with mctpye interpolate'
      sys.exit()
  else:
    pot = potpars
elif parmw == 'saddlepoint':
  potpars = np.random.rand(bins, paratyps)
  Params = np.genfromtxt(parname, skip_header=1)
  atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)
  if np.shape(Params)[-1] != atyps:
    print 'make sure that the parameter file has size ',atyps,'x',atyps
    sys.exit()
  ind = 0
  for n in range(atyps):
    for nn in range(n, atyps):
      potpars[0,ind] = Params[n,nn]
      potpars[1,ind] = Params[n+atyps,nn]
      potpars[-1,ind] = Params[n-atyps,nn]
      ind += 1
  if mctype == 'interpolate':
    pot = fgenlib.fgen(stepsize, start, rcutoff, potsteps, fgentype, potpars, power1, power2, paratyps, -1)
  else: 
    pot = potpars  
elif parmw == 'ABC':
  potpars = np.random.rand(bins, paratyps)
  Params = np.genfromtxt(parname, skip_header=1)
  rbins = np.shape(Params)[0]/np.shape(Params[1])
  atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)
  if np.shape(Params)[-1] != atyps:
    print 'make sure that the parameter file has size ',atyps,'x',atyps
    sys.exit()
  ind = 0
  for b in range(rbins):
    for n in range(atyps):
      for nn in range(n, atyps):
	potpars[b,ind] = Params[b*atyps+n,nn]
      ind += 1      
  pot = potpars
else:
  print 'no paremtertype given'
  sys.exit()

if dup == 'irmsd' or dup == 'lrmsd':
  weight[:,:]=1./rmsd[:,:]
  if '--cutoffweight' in sys.argv:
    evcut=1./evcut
  if '--maxweight' in sys.argv:
    maxcut = 1./maxcut
elif dup=='fnat' or dup=='capstars':
  weight[:,:]=rmsd[:,:]
elif dup == 'probabilities':
  targetvalues=np.genfromtxt(valuate)	#weights for targetfunction
  for i in range(len(targetvalues)-1):
    rmin=targetvalues[i,0]
    rmax=targetvalues[i+1,0]
    maske=ma.masked_outside(rmsd,rmin,rmax+0.00001)
    coor=ma.nonzero(maske)
    for j in range(len(coor[0])):
      weight[coor[0][j],coor[1][j]]=targetvalues[i,prorow]		#targetvalues are in the file 1-3 prob, 4 weight

if '--cutoffweight' in sys.argv:
    weight=ma.masked_less(weight,evcut).filled(0.)

if '--maxweight' in sys.argv:
  weight = ma.masked_greater(weight, maxcut).filled(maxcut)

#test if all the complexes have at least one good structure to train and evaluate on
norm=np.sum(weight, axis=1)		#sums all weight up to norm the evaluation for each complex
masknorm=ma.array(cmplx,mask=ma.getmask(ma.masked_not_equal(norm,0)))
masknorm=ma.compressed(masknorm)
if len(masknorm)!=0:
  print 'ATTENTION!!!: '+str(len(masknorm))+' complexes do not have any good structure to train on'
  print masknorm
  sys.exit()

reweight=np.zeros((numcpx,structures))
#choice of targetfunction
if target=='refine':
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=np.sum(reweight[:,:bestof], axis=1)/norm
      ausgabe=evaluo[:]
      return evaluo,ausgabe
if target=='rank':
    position=np.arange(bestof,0,-1)
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=np.sum(reweight[:,:bestof]*position,axis=1)/norm
      ausgabe=np.sum(reweight[:,:bestof], axis=1)/norm
      return evaluo,ausgabe
if target=='simple':			#usually use another valuate function
    position=np.arange(1,structures+1)
    def targfunc(eni,weighti):
      #evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)
      for i in range(numcpx):
        reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo = -(ma.notmasked_edges(ma.masked_equal(reweight,0),axis=1)[0][1]+1.)
      ausgabe=evaluo[:]
      return evaluo,ausgabe
if target=='positionlinear':
    position=np.arange(structures-1,-1,-1)
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=np.sum(reweight[:,:]*position,axis=1)/norm
      ausgabe=evaluo[:]
      return evaluo,ausgabe
if target=='refine-positionlinear':
    position=np.arange(structures-1,-1,-1)
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=(np.sum(reweight[:,:]*position,axis=1)/norm)*(np.sum(reweight[:,:bestof], axis=1)/norm)
      ausgabe=evaluo[:]
      return evaluo,ausgabe
if target=='positionquadratic':
    position=(0.01*(np.arange(structures-1,-1,-1)))**2
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=np.sum(reweight[:,:]*position,axis=1)/norm
      ausgabe=evaluo[:]
      return evaluo,ausgabe
if target=='positiontop':
    position=np.append(np.ones(bestof)*(structures-bestof),np.arange(structures-bestof,0,-1))    
    def targfunc(eni,weighti):
      evaluo=np.zeros(numcpx,dtype=np.float32)
      aus=np.argsort(eni, axis=-1)	#sorts position of energy after engery for every complex
      for i in range(numcpx):
	  reweight[i]=weighti[i][aus[i,:]]	#sorts the weights after the new positions of the energy
      evaluo=np.sum(reweight[:,:]*position,axis=1)/norm
      ausgabe=evaluo[:]
      return evaluo,ausgabe

#Monte-Carlo steps
steps=int(ansteps)			#steps of Monte Carlo Annealing
stepsfloat=float(steps)

if '--converge' in sys.argv:
  conlen = int(sys.argv[sys.argv.index('--converge')+1])
  maximumscore,maxausgabe=targfunc(-weight,weight)
  maximumnorm=np.sum(maximumscore[:(numcpx/crossnum)*(crossnum-1)])
  def converge(changevec):
    result=np.sum(changevec)/maximumnorm
    return result
else:
  conlen = 2
  def converge(changevec):
    return 1.

T0=float(temp)				#T0 has to be chosen in comparison to number of structures and the aimfunction
Emaxx=float(numcpx)*0.1

#choice of annealing scheme
if annealing=='linear':
    def function(k):
	return (T0)-(T0/stepsfloat)*k		#linear annealing
elif annealing=='exponential':
    def function(k):
	return T0*(0.001)**(k/(stepsfloat))		#exponential annealing
elif annealing=='logarithmic':
    def function(k):
	return Emaxx/np.log(2+k)			#logarithmic annealing
elif annealing=='ziczac':				#zic-zac annealing
    def function(k): 
	return T0*(0.001)**(k/(stepsfloat))*np.cos(2*np.pi*k/(stepsfloat*0.3))**2+(T0/10.)*(0.001)**(k/(stepsfloat))
elif annealing=='wallclimber':
    def function(k): 
	return 1
elif annealing=='constant':
    def function(k):
      return T0

if annealing=='wallclimber':
  def acfunc():
    return 1.
else:
  def acfunc():
    return random.random()

#choice of stepsize
if stepping=='normal':
  if mctype=='interpolate':
    if fgentype == 0:
      def stepfunc(st,temp,par):
	return (1+9*par)*random.choice((st,-st))
    else:
      def stepfunc(st,temp,par):
	return (-9*(par-1)**2+10)*random.choice((st,-st))
  else:
    def stepfunc(st,temp):
      return random.choice((st,-st))
elif stepping=='adaptive':
  if mctype=='interpolate':
    if fgentype == 0:
      def stepfunc(st,temp,par):
	return st*(1+9*par)*random.choice((temp/T0,temp/T0))
    else:
      def stepfunc(st,temp,par):
	return st*(-9*(par-1)**2+10)*random.choice((temp/T0,temp/T0))
  else:
    def stepfunc(st,temp):
      return st*random.choice((temp/T0,-temp/T0))

#MC precalculations
en=np.zeros((numcpx,structures),dtype=np.float32)
enold=np.zeros((numcpx,structures),dtype=np.float32)
for j in range(numcpx):    
    counts = allcounts[j]
    for i in range(structures):
	energy=np.sum(counts[:,i,:]*pot)
	en[j,i]=energy
	enold[j,i]=energy

evaluold,ausgab=targfunc(en,weight)
ausgeben=ausgab[:]

crossfile=open('Annealing-'+output+'.txt','w')
for i in cmplx:
  crossfile.write("%7s" % i)
crossfile.write('\n'+str(ausgeben).strip('[]'))


#start of MC-loop
#random.seed(0)
k=0
ac = 0

convec=np.ones(conlen)
son=np.arange(1,conlen)


if mctype=='interpolate':	
    potparsnew=np.zeros(len(potpars))
    if constpar :
      def changefunc(potpar, time, f):
	result=max(min(potpar+time, parranges[f,1]),parranges[f,0])
	return result
    else:
      def changefunc(potpar, time, f):
	result=potpar+time
	return result      
    while True:
	k+=1
	if k==steps or converge(convec)<1e-6: break 
	T=function(k)
	a=random.randint(0, bins-2)
	b=random.randint(0, paratyps-1)
	times=stepfunc(parstepsize,T,a)			#stepsize depends on what you want to change
	potparsnew[:]=potpars[:,b]
	potparsnew[a]=changefunc(potparsnew[a],times,a)
	potn = fgenlib.fgen(stepsize, start, rcutoff, potsteps, fgentype, potparsnew ,power1, power2, paratyps, 1)
	potnew = (np.ones((structures, potsteps))*potn).T
	potold = (np.ones((structures, potsteps))*pot[:,b]).T
	counts = allcounts[:,:,:,b]
	en[:,:]=enold+np.sum(counts*(potnew-potold), axis=1)
	evalu,ausgab=targfunc(en,weight)
	ins=-np.sum(evaluold[:(numcpx/crossnum)*(crossnum-1)]-evalu[:(numcpx/crossnum)*(crossnum-1)])
	test=np.exp(ins/T)
	acceptance=acfunc()
	if acceptance<test:
	    potpars[a,b]=changefunc(potpars[a,b],times,a)
	    potpars[-1,b]=np.sign(np.sign(potpars[1,b])+0.01)*potpars[-1,b]
	    potpars[a,b]=abs(potpars[a,b])
	    pot[:,b] = potn[:]
	    enold[:,:]=en[:,:]
	    evaluold[:]=evalu[:]
	    ausgeben[:]=ausgab[:]
	    ac+=1
	    convec=np.append(convec[son],[ins])
	else:
	    convec=np.append(convec[son],[0.]) 
	crossfile.write('\n'+str(ausgeben).strip('[]'))
    crossfile.close()
    
    atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)
    paramatrix=np.zeros((atyps*bins,atyps))
    ind=0
    for n in range(atyps):
	for nn in range(n,atyps):
	    for nnn in range(bins):
		paramatrix[nnn*atyps+n,nn]=potpars[nnn,ind]
		paramatrix[nnn*atyps+nn,n]=potpars[nnn,ind]
	    ind+=1
    np.savetxt('MC-Paramatrix_'+output+'_parm.par', paramatrix, fmt="%6s", header = Header+' accepted steps: '+str(ac)+'/'+str(k))

elif mctype == 'keepsign' or mctype == 'normal':
  
  if constpar:
    if mctype == 'keepsign':
      def changefunc(val, sub, f):
	return abs(max(min(val + sub, parranges[f,1]),parranges[f,0]))
    elif mctype == 'normal':
      def changefunc(val, sub, f):
	return max(min(val + sub, parranges[f,1]),parranges[f,0])
  else:
    if mctype == 'keepsign':
      def changefunc(val, sub, f):
	return abs(val+sub)
    elif mctype == 'normal':
      def changefunc(val, sub, f):
	return val + sub
    
  if changing == 'normal':
    def addfunc(eno, potpar, allcount, tim, a, b):
      counts = allcount[:,a,:,b]
      return eno+counts*(changefunc(potpar[a,b],times,a)-potpar[a,b])
  elif changing == 'saddle':
    def addfunc(eno, potpar, allcount, tim, a, b):
      counts = allcount[:,:,:,b]
      potnew = np.zeros(np.shape(potpar[:bins,b]))
      potnew[:]=potpar[:bins,b]
      print potnew, tim
      potnew[a]=changefunc(potpar[a,b],tim,a)
      print potnew
      return eno+counts[:,0,:]*((potnew[1]*potnew[0]**power1)-(potpars[1,b]*potpars[0,b]**power1))+counts[:,1,:]*((potnew[1]*potnew[0]**power2)-(potpars[1,b]*potpars[0,b]**power2))
  
  while True:
      k+=1
      if k==steps or converge(convec)<1e-6: break 
      T=function(k)
      a=random.randint(0,bins-1)
      b=random.randint(0,paratyps-1)
      times=stepfunc(parstepsize,T)
      en[:,:] = addfunc(enold, pot, allcounts, times, a, b) 
      evalu,ausgab=targfunc(en,weight)#score is the sum of the top 100 structures
      ins=-np.sum(evaluold[:(numcpx/crossnum)*(crossnum-1)]-evalu[:(numcpx/crossnum)*(crossnum-1)])
      test=np.exp(ins/T)
      acceptance=acfunc()
      if acceptance<test:
	  pot[a,b] = changefunc(pot[a,b],times,a)  ###potpars ????
	  enold[:,:]=en
	  evaluold[:]=evalu
	  ausgeben[:]=ausgab
	  ac+=1
	  convec=np.append(convec[son],[ins])
      else:
	 convec=np.append(convec[son],[0.]) 
      crossfile.write('\n'+str(ausgeben).strip('[]'))
  
  crossfile.close()
  
  if outform == 'matrix':
    outpar = np.zeros((bins, paratyps))
    
    if changing == 'saddle':
      outpar[0] = pot[1]*pot[0]**power1
      outpar[1] = pot[1]*pot[0]**power2
    else:
      outpar[:] = pot[:]
    
    if mctype == 'keepsign':   ##new convention
      for nnn in range(bins):
	if (nnn+1)%2 == 0:
	  outpar[nnn]=-outpar[nnn]
	  
    atyps = int(np.sqrt(0.25+2.*paratyps)-0.5)
    paramatrix=np.zeros((atyps*bins,atyps))
    ind=0
    for n in range(atyps):
      for nn in range(n,atyps):
	for nnn in range(bins):
	  paramatrix[nnn*atyps+n,nn]=outpar[nnn,ind]
	  paramatrix[nnn*atyps+nn,n]=outpar[nnn,ind]
	ind+=1
    np.savetxt('MC-Paramatrix_'+output+'.par', paramatrix, fmt="%6s", header = Header+' accepted steps: '+str(ac)+'/'+str(k))

    if changing == 'saddle':
      paramatrix = np.ones((3*atyps, atyps))
      ind=0
      for n in range(atyps):
	for nn in range(n,atyps):
	  for nnn in range(bins):
	    paramatrix[nnn*atyps+n,nn]=pot[nnn,ind]
	    paramatrix[nnn*atyps+nn,n]=pot[nnn,ind]
	  ind+=1
      np.savetxt('MC-Paramatrix_'+output+'_parm.par', paramatrix, fmt="%6s", header = Header+' accepted steps: '+str(ac)+'/'+str(k))
  
  else:
    np.savetxt('MC-Parameter_'+output+'.par', pot, fmt="%6s", header = Header+' accepted steps: '+str(ac)+'/'+str(k))
      
