import sys, os
import numpy as np
import numpy.ma as ma


if '--supplement' in sys.argv:
  sup = True
  supi = 2
else:
  sup = False
  supi = 1

if '--grids' in sys.argv:
  fullgrid=sys.argv[sys.argv.index('--grids')+1]
  outgridname=os.path.splitext(fullgrid)[0]
  fullgrids=np.load(fullgrid,mmap_mode='r')
  end = np.shape(fullgrids)[-2]
  if sup:
    trueint=sys.argv[sys.argv.index('--grids')+2]
    trueinterface=np.load(trueint)
else:
  print 'please insert grids'
  sys.exit()

if '--sorton' in sys.argv:
  sortons = sys.argv[sys.argv.index('--sorton')+1]
  cutmin=-float(sys.argv[sys.argv.index('--sorton')+2])
  sorton = -np.genfromtxt(sortons)[:end]
  outnamesorton = os.path.splitext(sortons)[0]
  sortname = os.path.splitext(sortons)[1][1:]
  if sup:
    supsortons = sys.argv[sys.argv.index('--sorton')+3]
    supsorton = -np.genfromtxt(supsortons)
  if os.path.splitext(sortons)[1][1:6] == 'lrmsd' or os.path.splitext(sortons)[1][1:6] == 'irmsd':
    cutmin = -cutmin
    sorton = -sorton
    if sup:
      supsorton = -supsorton
else:	
  print 'please give features to sort on'
  sys.exit()

if '--extras' in sys.argv:
  numextra = int(sys.argv[sys.argv.index('--extras')+1])
  nameextra =[]
  extras = []
  for i in range(numextra):
    nameextra.append(sys.argv[sys.argv.index('--extras')+2+i*supi])
    if nameextra[i].split('.')[-1][1:5]=='rmsd':
      extras.append(np.genfromtxt(sys.argv[sys.argv.index('--extras')+2+i*supi])[:end,-1])
      if sup:
	extras.append(np.genfromtxt(sys.argv[sys.argv.index('--extras')+2+i*supi+1])[:end,-1])      
    else:
      extras.append(np.genfromtxt(sys.argv[sys.argv.index('--extras')+2+i*supi])[:end])
      if sup:
	extras.append(np.genfromtxt(sys.argv[sys.argv.index('--extras')+2+i*supi+1])[:end])

if '--deleteclashes' in sys.argv:
  clashes = 'noclash'
  fullenergies=open(sys.argv[sys.argv.index('--deleteclashes')+1],'r')
  lines=fullenergies.readlines()
  fullvan=[]
  for i,line in enumerate(lines):
    if line[3:9]=='Energy':
      if lines[i+1][9:15]!='******':
	fullvan.append(float(lines[i+1][8:15]))
      else:
	fullvan.append(100000)
  fullvan=np.array(fullvan)[:end]
  if sup:
    trueenergies=open(sys.argv[sys.argv.index('--deleteclashes')+2],'r')
    lines=trueenergies.readlines()
    truevan=[]
    for i,line in enumerate(lines):
      if line[3:9]=='Energy':
	if lines[i+1][9:15]!='******':
	  truevan.append(float(lines[i+1][8:15]))
	else:
	  truevan.append(100000)
    truevan=np.array(truevan)
else:       
  clashes ='withclash'
  
if '--numset' in sys.argv:
  full = int(sys.argv[sys.argv.index('--numset')+1])
else:
  full = 100000

##### delete good structures in the top x% or in the last -x% ####
delfraction = ''
if '--deletefraction' in sys.argv:
  fraction = float(sys.argv[sys.argv.index('--deletefraction')+1])
  if np.sign(fraction) == -1:
    delfraction = 'best'+str(int(100+fraction))
    beginfrac = int((100.+fraction)*end/100.)
    delfrac = ma.nonzero(ma.masked_greater(sorton[beginfrac:],cutmin))[0]+beginfrac
    fullgrids=np.delete(fullgrids,delfrac, axis=-2)
    sorton=np.delete(sorton,delfrac,axis=0)
    if '--extras' in sys.argv:
      extras = np.delete(extras,delfrac,axis=1)
  elif np.sign(fraction) == 1:
    delfraction = 'last'+str(int(100-fraction))
    beginfrac = int(fraction*end/100.)
    delfrac = ma.nonzero(ma.masked_greater(sorton[:beginfrac],cutmin))[0]
    fullgrids=np.delete(fullgrids,delfrac, axis=-2)
    sorton=np.delete(sorton,delfrac,axis=0)
    if '--extras' in sys.argv:
      extras = np.delete(extras,delfrac,axis=1)

  
###sort good structures in front of the grid

if '--sortorig' in sys.argv:
  resort = '-resorted'+sortname+str(cutmin)
  sorti=np.argsort(sorton[:])
  sorts=np.arange(len(sorti))
  length=ma.size(ma.nonzero(ma.masked_greater(sorton[sorti],cutmin)))
  if length==0:
    print 'no good structures found in orignal decoyset'
  sorting=np.append(sorti[:length], np.delete(sorts,sorti[:length]))

else:
  resort = ''
  sorting = np.arange(len(sorton))
  length = 0
  
outgrid=np.zeros(np.shape(fullgrids), dtype = np.float32)

#sort grid
if len(np.shape(fullgrids))==4:
    for i,rcut in enumerate(fullgrids):
	for j,bins in enumerate(rcut):
	    outgrid[i,j,:,:]=bins[sorting]
elif len(np.shape(fullgrids)) == 2:
  outgrid = fullgrids[sorting]
else:
    for i,bins in enumerate(fullgrids):
	outgrid[i,:,:]=bins[sorting]

outsorton = sorton[sorting]

if '--extras' in sys.argv:
  outextras = []
  for i in range(numextra):
    outextras.append(extras[i*supi][sorting])
  outextras = np.array(outextras)
    
if '--deleteclashes' in sys.argv:
  delfulvan = ma.nonzero(ma.masked_less(fullvan[sorting][:length],0.))
  if len(delfulvan)==length:
    print 'no good structures of original set resorted'
  outgrid=np.delete(outgrid,delfulvan, axis=-2)
  outsorton=np.delete(outsorton,delfulvan,axis=0)
  
  if '--extras' in sys.argv:
    outextras = np.delete(outextras,delfulvan,axis=1)

if len(np.shape(outgrid))==4:
  outgrid=outgrid[:,:,:full,:]
elif len(np.shape(outgrid))==3:
  outgrid=outgrid[:,:full,:]
elif len(np.shape(outgrid))==2:
  outgrid=outgrid[:full,:]
else: 
  print 'unknwon shape of grid'
  sys.exit()
  
outsorton=outsorton[:full]
if '--extras' in sys.argv:
  outextras= outextras[:,:full]

if sup == False:
    
    np.save(outgridname+resort+delfraction+'-'+str(full)+'-'+clashes+'.npy',outgrid)
    
    if sortname[:5] == 'lrmsd' or sortname[:5] == 'irmsd':
      np.savetxt(outnamesorton+resort+delfraction+'-'+str(full)+'-'+clashes+'.'+sortname,outsorton, fmt=["%s","%3.3f"])
    else:
      np.savetxt(outnamesorton+resort+delfraction+'-'+str(full)+'-'+clashes+'.'+sortname,-outsorton, fmt=["%3.3f"])
    if '--extras' in sys.argv:
      for i in range(numextra):
        if os.path.splitext(nameextra[i])[1][1:6] == 'irmsd' or os.path.splitext(nameextra[i])[1][1:6] == 'lrmsd':
	  np.savetxt(os.path.splitext(nameextra[i])[0]+resort+delfraction+'-'+str(full)+'-'+clashes+os.path.splitext(nameextra[i])[1],np.array(zip(np.arange(1,len(outextras[i])+1),outextras[i])), fmt=["%d","%3.3f"])
        else:
	  np.savetxt(os.path.splitext(nameextra[i])[0]+resort+delfraction+'-'+str(full)+'-'+clashes+os.path.splitext(nameextra[i])[1],outextras[i], fmt=["%3.3f"])
		
###add good structures from constrained grid to original grid

elif sup == True:
    maske=ma.nonzero(ma.masked_greater(supsorton,cutmin))[0]
    
    length2=len(maske)
    if length2==0:
      print 'no good structures found to be supplemented'
    
    if len(np.shape(trueinterface))==4:
	supplegrid=np.zeros((len(trueinterface),len(trueinterface[0]),len(maske), len(trueinterface[0,0,0])), dtype = np.float32)
	for i,rcut in enumerate(trueinterface):
	    for j,bins in enumerate(rcut):
		supplegrid[i,j,:,:]=bins[maske]
    elif len(np.shape(trueinterface)) == 2:
      supplegrid = trueinterface[maske]
    else:
	supplegrid=np.zeros((len(trueinterface),len(maske), len(trueinterface[0,0])))
	for i,bins in enumerate(trueinterface):
	    supplegrid[i,:,:]=bins[maske]

    supplesorton = supsorton[maske]
    
    if '--extras' in sys.argv:
      suppleextra = []
      for i in range(numextra):
	suppleextra.append(extras[i*supi+1][maske])
      suppleextra = np.array(suppleextra)
    
    
    if '--deleteclashes' in sys.argv:
      supplevan=truevan[maske]
      delsupvan=ma.nonzero(ma.masked_less(supplevan[:length2],0.))
      if len(delsupvan)==length2:
	print 'no good supplemented structures could be added'
      supplegrid=np.delete(supplegrid,delsupvan, axis=-2)
      supplesorton=np.delete(supplesorton,delsupvan,axis=0)
      suppleextra=np.delete(suppleextra,delsupvan,axis=1)
    
    ausgrid=np.append(supplegrid,outgrid,axis=-2)
    ausgrid = ausgrid.astype(np.float32)
    if len(np.shape(ausgrid))==4:
      ausgrid=ausgrid[:,:,:full,:]
    elif len(np.shape(ausgrid))==3:
      ausgrid=ausgrid[:,:full,:]
    elif len(np.shape(ausgrid))==2:
      ausgrid=ausgrid[:full,:]
    else: 
      print 'unknwon shape of grid'
      sys.exit()
    
    aussorton=np.append(supplesorton,outsorton,axis=0)[:full]
    if '--extras' in sys.argv:
      ausextra = np.append(suppleextra, outextras, axis = 1)[:,:full]
         
    np.save(outgridname+resort+'-'+str(full)+'-'+clashes+'-sup'+sortname+str(cutmin)+'.npy',ausgrid)
    if sortname[:5] == 'lrmsd' or sortname[:5] == 'irmsd':
      np.savetxt(outnamesorton+resort+'-'+str(full)+'-'+clashes+'-sup'+sortname+str(cutmin)+'.'+sortname,aussorton, fmt=["%s","%3.3f"])
    else:
      np.savetxt(outnamesorton+resort+'-'+str(full)+'-'+clashes+'-sup'+sortname+str(cutmin)+'.'+sortname,-aussorton, fmt=["%3.3f"])
    for i in range(numextra):
      if os.path.splitext(nameextra[i])[1][1:6] == 'irmsd' or os.path.splitext(nameextra[i])[1][1:6] == 'lrmsd':
	np.savetxt(os.path.splitext(nameextra[i])[0]+resort+'-'+str(full)+'-'+clashes+'-sup'+sortname+str(cutmin)+os.path.splitext(nameextra[i])[1],np.array(zip(np.arange(1,len(ausextra[i])+1), ausextra[i])), fmt=["%d","%3.3f"])
      else:
	np.savetxt(os.path.splitext(nameextra[i])[0]+resort+'-'+str(full)+'-'+clashes+'-sup'+sortname+str(cutmin)+os.path.splitext(nameextra[i])[1],ausextra[i], fmt=["%3.3f"])
 
