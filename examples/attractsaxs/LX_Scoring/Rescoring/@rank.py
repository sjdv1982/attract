# @rank.py
# python --input input.dat --proteinmodel (tobi,opls,attract) receptor.pdb ligang.pdb --vdwparams vdwparameter.txt [for opls and attract insert: cutoff(float,'None'); for tobi insert bins: --bins #bins edges_i] --gbscore receptor-aa.pdb ligand-aa.pdb --solvation solvationparameter.txt --combination (linear,probability,desicionfunction) combineparameter.txt(.pkl)
#optional --rerank to resort output.dat, --output to name outputfile (default: input-rescored.dat or input-resorted.dat)
# for nonlinear combination (probability, decisionfunction) --preprocessing could be important
#for output file, scores are ordered: vdwscore| gbscore| solvationscore
import numpy as np
import sys, os
import numpy.ma as ma
import collectlibpy as collectlib
import scorelib as sl

if len(sys.argv)<2:
  print 'use as follows:','\n', 'python --input input.dat --proteinmodel (tobi,opls,attract) receptor.pdb ligang.pdb --vdwparams vdwparameter.txt [for opls and attract insert: cutoff(float,"None"); for tobi insert bins: --bins #bins edges_i] --gbscore receptor-aa.pdb ligand-aa.pdb --solvation solvationparameter.txt --combination (linear,probability,desicionfunction) combineparameter.txt(.pkl)'
  sys.exit()
  
if '--input' in sys.argv:
  datfile=sys.argv[sys.argv.index('--input')+1]
  gobj=open(datfile)
  glines=gobj.readlines()
  structures=0
  prot=0
  for gline in glines:
    if gline[:6]=='#pivot':
      prot+=1
    if gline.count('#')==0:
      structures+=1
  structures=structures/2
  prot=max(prot,2)
  ensnum=np.zeros((structures,prot),dtype=int)
else:
  print '#### Insert input-files as follows: --input out.dat'
  sys.exit()

if '--proteinmodel' in sys.argv:
  modeltype=sys.argv[sys.argv.index('--proteinmodel')+1]
  receptor=sys.argv[sys.argv.index('--proteinmodel')+2]
  ligand=sys.argv[sys.argv.index('--proteinmodel')+3]
  collectinsert=[datfile,receptor,ligand]
else:
  print '#### Insert receptor and ligand models as follows: --proteinmodel modeltype(tobi,opls,attract) receptor.pdb ligand.pdb'
  sys.exit()



ensfiles=[]
gbensfiles = []
modefile = None
name =None
modeon= False
enson=False
anr = 0
while 1:
    anr += 1
        
    if anr > len(sys.argv)-1: break  
    arg = sys.argv[anr]

    if anr <= len(sys.argv)-3 and arg == "--ens":
      enson = True
      ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
      sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
      anr -= 3
      continue

    if anr <= len(sys.argv)-2 and arg == "--modes":
      modeon = True
      modefile = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue

    if anr <= len(sys.argv)-2 and arg == "--name":
      modeon = True
      name = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue

#if namefile: collectinsert+=['--name', namefile]
    
if modefile: collectinsert+=['--modes', modefile]

for nr, ensfile in ensfiles:
  collectinsert += ["--ens", nr, ensfile]
  cont=0
  ind=0
  for gline in glines:
    if gline.count('#')==0:
      cont+=1
      if cont%prot==float(nr)%prot:
	ensnum[ind,int(nr)-1]=int(gline.split()[0])
	ind+=1

gobj.close()

collectlib.collect_init(collectinsert)

if modeltype=='tobi':
  natyps=19
elif modeltype=='opls':
  natyps=13
elif modeltype=='attract':
  natyps=32
elif modeltype=='gaa':
  natyps=27
elif modeltype=='any':
  natyps=int(sys.argv[sys.argv.index('--proteinmodel')+4])
else:
  print 'no model for atomtypes chosen! use tobi, opls, gaa, attract or "any"!'

#put rvdwsol somewhere else to make vdwscoring faster
Rrvdwsol=[]
atyper=[]
fobj=open(receptor, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
    atyper.append(int(fline[57:59]))
fobj.close()
atyper=np.array(atyper)

Rlvdwsol=[]
atypel=[]
fobj=open(ligand, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
    atypel.append(int(fline[57:59]))
fobj.close()

atypel=np.array(atypel)

rlen = len(atyper)
llen = len(atypel)

if modeltype=='opls':
  mp = np.arange(0,83)
  replace={30:2, 31:3, 65:4, 66:5, 67:6, 68:7, 69:8, 70:9, 71:10, 80:11, 81:12, 82:13}
  mp[replace.keys()] = replace.values()
  atypel = mp[atypel]
  atyper = mp[atyper]
elif modeltype=='tobi':
  mp = np.arange(0,33)
  replace={32:19}
  mp[replace.keys()] = replace.values()
  atypel = mp[atypel]
  atyper = mp[atyper]
elif modeltype == 'attract':
  mp = np.arange(0,100)
  replace={99:32}
  mp[replace.keys()] = replace.values()
  atypel = mp[atypel]
  atyper = mp[atyper]


atypestruc=np.append(atyper,atypel)


if '--shift' in sys.argv:
  shift = float(sys.argv[sys.argv.index('--shift')+1])
else:
  shift = 0. 

combnumber=np.ones(5,dtype=np.int)
numscores=0
stepon=False
vdwon=False
solon=False
elon=False
bsaon=False

if '--vdwpotential' in sys.argv:
    vdwon=True
    vdwparams=np.genfromtxt(sys.argv[sys.argv.index('--vdwpotential')+1],skip_header=1)
    testmodel=len(vdwparams[0])
    
    if '--vdwcutoff' in sys.argv:
      vdwcutoff=float(sys.argv[sys.argv.index('--vdwcutoff')+1])
    else:
      vdwcutoff = 100.   
      
    if '--vdwfunctiontype' in sys.argv:
      vdwfunctiontype=sys.argv[sys.argv.index('--vdwfunctiontype')+1]
    else:
      print 'give type of vdw potential'
      sys.exit()
      
    sig = np.zeros((natyps,natyps))
    eps = np.zeros((natyps,natyps))
    ivor = np.ones((natyps,natyps))
    hexp = np.zeros((natyps,natyps))
    ratshift = np.zeros((natyps,natyps)) 
    
    ivor[:testmodel,:testmodel]=vdwparams[-1*testmodel:,:]
    sig[:testmodel,:testmodel]=vdwparams[:testmodel]
    eps[:testmodel,:testmodel]=vdwparams[testmodel:2*testmodel]
    
    exp = False
    atshift = False

    if vdwfunctiontype=='attractpartial':
      power1 = 8.
      power2 = 6.
      exp = True
	
    elif vdwfunctiontype=='attractlong':
      power1 = 8.
      power2 = 6.
      atshift = True
	
    elif vdwfunctiontype=='attract':
      power1 = 8.
      power2 = 6.
	
    elif vdwfunctiontype=='opls':
      power1 = 12.
      power2 = 6.
	
    elif vdwfunctiontype == 'free':
      power1 = float(sys.argv[sys.argv.index('--functiontype')+2])
      power2 = float(sys.argv[sys.argv.index('--functiontype')+3])

      if len(sys.argv)-sys.argv.index('--functiontype') > 4:
	if sys.argv[sys.argv.index('--functiontype')+4] == 'exp':
	  exp = True
	elif sys.argv[sys.argv.index('--functiontype')+4] == 'atshift':
	  atshift = True
      if len(sys.argv)-sys.argv.index('--functiontype') > 5:
	if sys.argv[sys.argv.index('--functiontype')+5] == 'exp':
	  exp = True
	elif sys.argv[sys.argv.index('--functiontype')+5] == 'atshift':
	  atshift = True
    else:
      print ' type of function not understood, please use one of the given forms or define your own function'
      sys.exit()

    if exp == True and atshift == False:
      hexp[:parcount, :parcount] = Params[parcount*2:parcount*3]
    elif exp == False and atshift == True:
      ratshift[:parcount, :parcount] = Params[parcount*2:parcount*3]
    elif exp == True and atshift == True:
      ratshift[:parcount, :parcount] = Params[parcount*2:parcount*3]
      hexp[:parcount, :parcount] = Params[parcount*3:parcount*4]

    def vdwscore(reccoor, ligcoor):
      result = sl.score(reccoor, ligcoor, shift, vdwcutoff, power1, power2, eps, sig, ivor, hexp, ratshift, atyper, atypel) 
      return result

#  if testmodel!=natyps:
#    print ' Vdw parameter do not match with the chosen protein model!'
#    sys.exit()
    combnumber[0]=numscores
    numscores+=1
    


if '--steppotential' in sys.argv:
    stepon=True
    combnumber[1]=numscores
    numscores+=1
    
    stepparams=np.genfromtxt(sys.argv[sys.argv.index('--steppotential')+1],skip_header=1)
    testmodel=len(stepparams[0])
    
    if '--bins' in sys.argv:
      bins=int(sys.argv[sys.argv.index('--bins')+1])
      cutoff=np.zeros(bins)
      steppar=np.zeros((bins, natyps, natyps))
      for i in range(bins):
	cutoff[i]=float(sys.argv[sys.argv.index('--bins')+2+i])**2
	steppar[i]=stepparams[i*testmodel:(i+1)*testmodel]
	#cutoff quadrieren

      def stepscoring(reccoor, ligcoor):
	en=sl.stepscore(cutoff, reccoor, ligcoor, atyper, atypel, steppar)  
	return en
    else:
      print 'if modeltype steps, one has to give the number of steps and the values of the edges wiht --bins #steps edge1 ...'
      sys.exit()
      


if '--solvation' in sys.argv:    
  solon=True
  solparams=np.genfromtxt(sys.argv[sys.argv.index('--solvation')+1],skip_header=1)
  combnumber[2]=numscores
  numscores+=1
  if len(solparams)!=natyps:
    print 'solvation parameter do not match with the chosen protein model!'
    sys.exit()
  import asalib
  
  #put rvdwsol somewhere else to make vdwscoring faster
  Rrvdwsol=[]
  fobj=open(receptor, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      if fline[13]=='N':
	Rrvdwsol.append(1.6)
      elif fline[13]=='C':
	Rrvdwsol.append(1.7)
      elif fline[13]=='O':
	Rrvdwsol.append(1.5)
      elif fline[13]=='S':
	Rrvdwsol.append(2.0)
      elif fline[13]=='H' or fline[12]=='H':
	Rrvdwsol.append(0.)
      elif fline[13]=='P':
	Rrvdwsol.append(2.0)
      else:
	print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', receptor
	sys.exit()
  fobj.close()
  #test=np.array(test)
  Rrvdwsol=np.array(Rrvdwsol)
  atyper=np.array(atyper)
  
  Rlvdwsol=[]
  fobj=open(ligand, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      if fline[13]=='N':
	Rlvdwsol.append(1.6)
      elif fline[13]=='C':
	Rlvdwsol.append(1.7)
      elif fline[13]=='O':
	Rlvdwsol.append(1.5)
      elif fline[13]=='S':
	Rlvdwsol.append(2.0)
      elif fline[13]=='H' or fline[12]=='H':
	Rlvdwsol.append(0.)
      elif fline[13]=='P':
	Rlvdwsol.append(2.0)
      else:
	print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', ligand
	sys.exit()
  fobj.close()
  Rlvdwsol=np.array(Rlvdwsol)
  atypel=np.array(atypel)
  
  if '--watersize' in sys.argv:
    wrad = float(sys.argv[sys.argv.index('--watersize')+1])
    wout = 'wrad'+str(wrad)
  else:
    wrad = 1.4
    wout = ''
  
  cgname = ''
  
  if '--coarse_grained' in sys.argv:
    cgname = '-cg-radii'
    mp = np.arange(0.,33.)
    replace={1:2.0, 2:1.9, 3:1.95, 4:1.9, 5:1.9, 6:1.9, 7:2.0, 8:2.0, 9:1.9, 10:2.0, 11:1.9, 12:2.0, 13:1.9, 14:2.2, 15:2.2, 16:2.0, 17:1.9, 18:2.0, 19:2.0, 20:2.2, 21:2.2, 22:1.9, 23:1.9, 24:1.9, 25:2.0, 26:2.2, 27:2.2, 28:2.2, 29:2.0, 30:1.6, 31:1.5, 32:1.7}
    mp[replace.keys()] = replace.values()
    Rrvdwsol = mp[atyper]
    Rlvdwsol = mp[atypel]
  
  Rvdwstrucsol=np.append(Rrvdwsol,Rlvdwsol)  
  enrsol=[]
  enlsol=[]
  coor = np.array(collectlib.collect_all_coor())
  coor_r=coor[:rlen]
  coor_l=coor[rlen:]
  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdwsol, atyper, natyps,wrad)
  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdwsol, atypel, natyps,wrad)
  enrsol.append(np.sum(asar*solparams))
  enlsol.append(np.sum(asal*solparams))
  if enson:
    for model in ensfiles:
      gobj=open(model[1],'r')
      glines=gobj.readlines()
      nummod=len(glines)
      if model[0]=='1':
	for remod in glines:
	  fobj=open(remod.strip())
	  coor=[]
	  for flines in fobj:
	    if flines[:4]=='ATOM':
	      coor.append([float(flines[30:38]),float(flines[38:46]),float(flines[46:54])])
	  coor=np.array(coor)
	  enrsol.append(np.sum(asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rrvdwsol, atyper, natyps)*solparams))
      if model[0]=='2':
	for remod in glines:
	  fobj=open(remod.strip())
	  coor=[]
	  for flines in fobj:
	    if flines[:4]=='ATOM':
	      coor.append([float(flines[30:38]),float(flines[38:46]),float(flines[46:54])])
	  coor=np.array(coor)
	  enlsol.append(np.sum(asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rlvdwsol, atypel, natyps)*solparams))
  enrsol.append(0.0)
  enlsol.append(0.0)
  
    
  if modeon:
    def solvatescore(rcoor, lcoor, coors):
      asar=asalib.asaatom(rcoor[:,0],rcoor[:,1],rcoor[:,2], Rrvdwsol, atyper, natyps,wrad)
      asal=asalib.asaatom(lcoor[:,0],lcoor[:,1],lcoor[:,2], Rlvdwsol, atypel, natyps,wrad)
      enrsolv=np.sum(asar*solparams)
      enlsolv=np.sum(asar*solparams)
      asastruc=asalib.asaatom(coors[:,0],coors[:,1],coors[:,2], Rvdwstrucsol, atypestruc, natyps,wrad)
      en = enrsolv+enlsolv-np.sum(asastruc*solparams)
      return en
  else:
    def solvatescore(rcoor, lcoor, coors):
      asastruc=asalib.asaatom(coors[:,0],coors[:,1],coors[:,2], Rvdwstrucsol, atypestruc, natyps,wrad)
      en = enrsol[ensnum[nstruc,0]]+enlsol[ensnum[nstruc,1]]-np.sum(asastruc*solparams)
      return en



if '--electrostatics' in sys.argv:
  elon=True
  charger=[]
  fobj=open(receptor, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      charger.append(float(fline[60:67]))
  fobj.close()
  #test=np.array(test)
  charger=np.array(charger)

  chargel=[]
  fobj=open(ligand, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      chargel.append(float(fline[60:67]))
  fobj.close()
  chargel=np.array(chargel)

  rlenmod = len(charger)
  llenmod = len(chargel)

  chargestruc=np.append(charger,chargel)  
  def elscore(rcoor, lcoor):
    en=sl.elec(coor_r, coor_l, charger, chargel, shift)
    return en
  combnumber[3]=numscores
  numscores+=1


if '--buriedsa' in sys.argv:
  bsaon=True
  combnumber[4]=numscores
  numscores+=1
  import asalib
  
  #put rvdwsol somewhere else to make vdwscoring faster
  Rrvdwsol=[]
  fobj=open(receptor, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      if fline[13]=='N':
	Rrvdwsol.append(1.6)
      elif fline[13]=='C':
	Rrvdwsol.append(1.7)
      elif fline[13]=='O':
	Rrvdwsol.append(1.5)
      elif fline[13]=='S':
	Rrvdwsol.append(2.0)
      elif fline[13]=='H' or fline[12]=='H':
	Rrvdwsol.append(0.)
      elif fline[13]=='P':
	Rrvdwsol.append(2.0)
      else:
	print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', receptor
	sys.exit()
  fobj.close()
  #test=np.array(test)
  Rrvdwsol=np.array(Rrvdwsol)
  
  Rlvdwsol=[]
  fobj=open(ligand, 'r')
  for i,fline in enumerate(fobj):
    if fline[:4]=='ATOM':
      if fline[13]=='N':
	Rlvdwsol.append(1.6)
      elif fline[13]=='C':
	Rlvdwsol.append(1.7)
      elif fline[13]=='O':
	Rlvdwsol.append(1.5)
      elif fline[13]=='S':
	Rlvdwsol.append(2.0)
      elif fline[13]=='H' or fline[12]=='H':
	Rlvdwsol.append(0.)
      elif fline[13]=='P':
	Rlvdwsol.append(2.0)
      else:
	print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', ligand
	sys.exit()
  fobj.close()
  Rlvdwsol=np.array(Rlvdwsol)
  
  if '--watersize' in sys.argv:
    wrad = float(sys.argv[sys.argv.index('--watersize')+1])
    wout = 'wrad'+str(wrad)
  else:
    wrad = 1.4
    wout = ''
  
  cgname = ''
  
  if '--coarse_grained' in sys.argv:
    cgname = '-cg-radii'
    mp = np.arange(0.,33.)
    replace={1:2.0, 2:1.9, 3:1.95, 4:1.9, 5:1.9, 6:1.9, 7:2.0, 8:2.0, 9:1.9, 10:2.0, 11:1.9, 12:2.0, 13:1.9, 14:2.2, 15:2.2, 16:2.0, 17:1.9, 18:2.0, 19:2.0, 20:2.2, 21:2.2, 22:1.9, 23:1.9, 24:1.9, 25:2.0, 26:2.2, 27:2.2, 28:2.2, 29:2.0, 30:1.6, 31:1.5, 32:1.7}
    mp[replace.keys()] = replace.values()
    Rrvdwsol = mp[atyper]
    Rlvdwsol = mp[atypel]
  
  Rvdwstrucsol=np.append(Rrvdwsol,Rlvdwsol)  
  
  
  enrbsa=[]
  enlbsa=[]
  coor = np.array(collectlib.collect_all_coor())
  coor_r=coor[:rlen]
  coor_l=coor[rlen:]
  enrbsa.append(asalib.asa(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdwsol, wrad))
  enlbsa.append(asalib.asa(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdwsol,wrad))
  if enson:
    for model in ensfiles:
      gobj=open(model[1],'r')
      glines=gobj.readlines()
      nummod=len(glines)
      if model[0]=='1':
	for remod in glines:
	  fobj=open(remod.strip())
	  coor=[]
	  for flines in fobj:
	    if flines[:4]=='ATOM':
	      coor.append([float(flines[30:38]),float(flines[38:46]),float(flines[46:54])])
	  coor=np.array(coor)
	  enrbsa.append(asalib.asa(coor[:,0],coor[:,1],coor[:,2], Rrvdwsol, wrad))
      if model[0]=='2':
	for remod in glines:
	  fobj=open(remod.strip())
	  coor=[]
	  for flines in fobj:
	    if flines[:4]=='ATOM':
	      coor.append([float(flines[30:38]),float(flines[38:46]),float(flines[46:54])])
	  coor=np.array(coor)
	  enlbsa.append(asalib.asa(coor[:,0],coor[:,1],coor[:,2], Rlvdwsol, wrad))
  enrsol.append(0.0)
  enlsol.append(0.0)
  
    
  if modeon:
    def bsa(rcoor, lcoor, coors):
      enrbs=asalib.asa(rcoor[:,0],rcoor[:,1],rcoor[:,2], Rrvdwsol, wrad)
      enlbs=asalib.asa(lcoor[:,0],lcoor[:,1],lcoor[:,2], Rlvdwsol,wrad)
      asastruc=asalib.asa(coors[:,0],coors[:,1],coors[:,2], Rvdwstrucsol, wrad)
      en = -(enrbs+enlbs-asastruc)
      return en
  else:
    def bsa(rcoor, lcoor, coors):
      asastruc=asalib.asa(coors[:,0],coors[:,1],coors[:,2], Rvdwstrucsol, wrad)
      en = -(enrbsa[ensnum[nstruc,0]]+enlbsa[ensnum[nstruc,1]]-asastruc)
      return en  


    
if numscores==0:
  print ' No parameter given: insert parameter with --vdwparams, --solvation or --gbscore'
  sys.exit()

combineparams=np.ones(5)  

if '--combination' in sys.argv:
  method=params=sys.argv[sys.argv.index('--combination')+1]
  if method=='linear':
    combparams=np.genfromtxt(sys.argv[sys.argv.index('--combination')+2],skip_header=1)
    if len(combparams)!=numscores:
      print ' number of parameter to combine: ',len(combparams)
      print ' number of parameter given: ', numscores
      sys.exit()
    else:
	combineparams=combparams[combnumber]
    def combinescores(inscores):
      endscores=inscores*combineparams
      endscore=np.sum(endscores,axis=1)
      return endscores,endscore
    
  elif method=='probability':
    import cPickle
    with open(combineparams, 'rb') as fid:
      wclf = cPickle.load(fid)
    if 'preprocessing' in sys.argv and (sys.argv[sys.argv.index('--preprocessing')+1]!='Standard' or sys.argv[sys.argv.index('--preprocessing')+1]!='MinMax'):
      print '#### no scaler chosen for preprocessing ####'
      sys.exit()
    def combinescores(inscores):
      endscore=wclf.predict_proba(inscores)
      return inscores, endscore
    
  elif method=='decisionfunction':
    import cPickle
    with open(combineparams, 'rb') as fid:
      wclf = cPickle.load(fid)
    if '--preprocessing' in sys.argv and (sys.argv[sys.argv.index('--preprocessing')+1]!='Standard' or sys.argv[sys.argv.index('--preprocessing')+1]!='MinMax'):
      print '#### no scaler chosen for preprocessing ####'
      sys.exit()
    def combinescores(inscores):
      endscore=wclf.decision_function(inscores)
      return inscores, endscore
  else:
    print 'no method for combination given! Insert combination parameter as follows: --combination method(linear,nonlinear) parameter.txt(pkl)'
    sys.exit()
else:
  def combinescores(inscores):
      endscores=inscores*combineparams
      endscore=np.sum(endscores,axis=1)
      return endscores,endscore

outscores=[]  
nstruc=0

while 1:
  scores=np.zeros(5)
  
  if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
  
  result = collectlib.collect_next()
  if result: break
  coor = np.array(collectlib.collect_all_coor())
  coor_r = coor[:rlen]
  coor_l = coor[rlen:]
  if vdwon:
    scores[0]=vdwscore(coor_r, coor_l)
  if stepon:
    scores[1]=stepscoring(coor_r, coor_l)
  if solon:
    scores[2]=solvatescore(coor_r, coor_l,coor)
  if elon:
    scores[3]=elscore(coor_r, coor_l)
  if bsaon:
    scores[4]=bsa(coor_r, coor_l,coor)
  outscores.append(scores)
  nstruc+=1

outscores=np.array(outscores)

    


if '--preprocessing' in sys.argv:
  from sklearn import preprocessing
  if sys.argv[sys.argv.index('--preprocessing')+1]=='Standard':
      scaler=preprocessing.StandardScaler().fit(outscores)
      outscores=scaler.transform(outscores)*100.
  elif sys.argv[sys.argv.index('--preprocessing')+1]=='MinMax':
      scaler=preprocessing.MinMaxScaler().fit(outscores)
      outscores=scaler.transform(outscores)*100.

#combine scores
outscores,finalscore = combinescores(outscores)


#print finalscore
if '--printout' in sys.argv:
  for fscore in finalscore:
    print fscore
    
if '--output' in sys.argv:
  outputname=sys.argv[sys.argv.index('--output')+1]
elif '--output' not in sys.argv and '--rescore' in sys.argv: 
  outputname=os.path.splitext(datfile)[0]+'-rescored.dat'
elif '--output' not in sys.argv and '--rerank' in sys.argv:
  outputname=os.path.splitext(datfile)[0]+'-resorted.dat'

obj = open(datfile,'r')
outlines=obj.readlines()
obj.close()

for outline in outlines:
  if outline[:6]=='#pivot':
    break
  else:
    outlines=outlines[1:]

    
if '--rerank' in sys.argv:
  outobj=open(outputname,'w')
  pivot=[]
  oldnum=[]
  seed=[]
  der=[]
  grad=[]
  degree=[]
  t=0
  for i,line in enumerate(outlines):
    
    if i<4:
      pivot.append(line)
    if i>=4 and line.count('#')==1:
      oldnum.append(line[1:].strip())
    if line[:8]=='### SEED':
      seed.append(line)
    if line[-12:]=='deredundant\n':
      der.append(line)
    if (line.count('#')==0) and (t%2==0):
      t+=1
      grad.append(line)
    elif (line.count('#')==0) and (t%2==1):
      t+=1
      degree.append(line)
      
  sort=np.argsort(finalscore)
  oldnum=np.array(oldnum,dtype=str)[sort]
  finalscore=finalscore[sort]
  outscores=outscores[sort]
  seed=np.array(seed,dtype=str)[sort]
  grad=np.array(grad,dtype=str)[sort]
  degree=np.array(degree,dtype=str)[sort]

  if '--numsubset' in sys.argv:
    nstruc=int(sys.argv[sys.argv.index('--numsubset')+1])

  for piv in pivot:
    outobj.write(piv)
  if len(der)!=0:
    der=np.array(der,dtype=str)[sort]
    for i in range(nstruc):
      outobj.write('#'+str(i+1)+'\n'+seed[i]+'##'+oldnum[i]+' => sort\n'+'## Energy: '
		  +"%.9f" % finalscore[i]+'\n'+'## '+"%12s" % str(round(outscores[i,0],3))
		  + "%12s" % str(round(outscores[i,1],3))+ "%12s" % str(round(outscores[i,2],3))
		  + "%12s" % str(round(outscores[i,3],3))+ "%12s" % str(round(outscores[i,4],3))+ "%12s" % '0.000'+'\n'+der[i]+grad[i]+degree[i])
  else:
    for i in range(nstruc):
      outobj.write('#'+str(i+1)+'\n'+seed[i]+'##'+oldnum[i]+' => sort\n'+'## Energy: '
		  +"%.9f" % finalscore[i]+'\n'+'## '+"%12s" % str(round(outscores[i,0],3))
		  + "%12s" % str(round(outscores[i,1],3))+ "%12s" % str(round(outscores[i,2],3))
		  + "%12s" % str(round(outscores[i,3],3))+ "%12s" % str(round(outscores[i,4],3))+ "%12s" % '0.000'+'\n'+grad[i]+degree[i])
  outobj.close()
elif '--rescore' in sys.argv:
  outobj=open(outputname,'w')
  j=0
  for i in range(len(outlines)):
    if outlines[i][:10]=='## Energy:':
      outobj.write('## Energy: '+"%6.9f" % finalscore[j]+'\n')
      outobj.write('## '+"%12s" % str(round(outscores[j,0],3))+ "%12s" % str(round(outscores[j,1],3))+ "%12s" % str(round(outscores[j,2],3))+ "%12s" % str(round(outscores[j,3],3))+ "%12s" % str(round(outscores[j,4],3))+ "%12s" % '0.000'+'\n')
      j+=1
    elif outlines[i-1][:10]!='## Energy:':
      outobj.write(outlines[i])
  outobj.close()


  
  
  
  
  
