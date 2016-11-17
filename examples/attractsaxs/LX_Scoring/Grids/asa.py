#collect.py
import numpy as np
import numpy.ma as ma
import sys,os
from scipy.spatial.distance import cdist
import collectlibpy as collectlib
import asalib

receptor=sys.argv[1]
ligand=sys.argv[2]
data=sys.argv[3]
surface=sys.argv[4]

if '--watersize' in sys.argv:
  wrad = float(sys.argv[sys.argv.index('--watersize')+1])
  wout = 'wrad'+str(wrad)
else:
  wrad = 1.4
  wout = ''

collectinsert=[data,receptor,ligand]

modefile=None
name=None
modeon=False
ensfiles=[]
anr = 0
while 1:
    anr += 1
        
    if anr > len(sys.argv)-1: break  
    arg = sys.argv[anr]

    if anr <= len(sys.argv)-3 and arg == "--ens":
      ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
      sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
      anr -= 3
      continue

    if anr <= len(sys.argv)-2 and arg == "--modes":
      modeon=True
      modefile = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue
    
    if anr <= len(sys.argv)-2 and arg == "--name":
      modeon=True
      name = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue

if modefile: collectinsert+=['--modes', modefile]
#if namefile: collectinsert+=['--name', namefile]

for nr, ensfile in ensfiles:
  collectinsert += ["--ens", nr, ensfile]

collectlib.collect_init(collectinsert)

atypr=[]
modatypr=[]
Rrvdw=[]
charger=[]
#test=[]
fobj=open(receptor, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
#    test.append([float(fline[30:38]),float(fline[38:46]),float(fline[46:54])])
    charger.append(float(fline[60:67]))
    modatypr.append(int(fline[57:59]))
    if fline[13]=='N':
      Rrvdw.append(1.6)
      atypr.append(2)
    elif fline[13]=='C':
      Rrvdw.append(1.7)
      atypr.append(1)
    elif fline[13]=='O':
      Rrvdw.append(1.5)
      atypr.append(3)
    elif fline[13]=='S':
      Rrvdw.append(2.0)
      atypr.append(4)
    elif fline[13]=='H' or fline[12]=='H':
      Rrvdw.append(0.)
      atypr.append(0)
    elif fline[13]=='P':
      Rrvdw.append(2.0)
      atypr.append(5)
    else:
      print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', receptor
      sys.exit()
fobj.close()
#test=np.array(test)
Rrvdw=np.array(Rrvdw)
charger=np.array(charger)
atypr=np.array(atypr)
modatypr=np.array(modatypr)

atypl=[]
modatypl=[]
Rlvdw=[]
chargel=[]
fobj=open(ligand, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
    chargel.append(float(fline[60:67]))
    modatypl.append(int(fline[57:59]))
    if fline[13]=='N':
      Rlvdw.append(1.6)      
      atypl.append(2)
    elif fline[13]=='C':
      Rlvdw.append(1.7)
      atypl.append(1)
    elif fline[13]=='O':
      Rlvdw.append(1.5)
      atypl.append(3)
    elif fline[13]=='S':
      Rlvdw.append(2.0)
      atypl.append(4)
    elif fline[13]=='H' or fline[12]=='H':
      Rlvdw.append(0.)
      atypl.append(0)
    elif fline[13]=='P':
      Rlvdw.append(2.0)
      atypl.append(5)
    else:
      print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', ligand
      sys.exit()
fobj.close()
Rlvdw=np.array(Rlvdw)
chargel=np.array(chargel)
atypl=np.array(atypl)
modatypl=np.array(modatypl)

Rvdwstruc=np.append(Rrvdw,Rlvdw)
chargestruc=np.append(charger,chargel)
atypstruc=np.append(atypr, atypl)
modatypstruc=np.append(modatypr, modatypl)

rlen = len(charger)
nstruc=0
#print 'receptor ', asar
#print 'ligand ', asal
if surface=='buriedsa':
    
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asa(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw,wrad)
      asal=asalib.asa(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
      
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon==True:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asa(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw,wrad)
	  asal=asalib.asa(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw,wrad)
	asa=asalib.asa(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc,wrad)
	print -(asar+asal-asa)
elif surface=='asa':
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	asa=asalib.asa(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc,wrad)
	print asa
elif surface=='atomsurface':
    natyps=np.amax(atypstruc)
    grid=[]
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, atypr, natyps,wrad)
      asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, atypl, natyps,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, atypr, natyps,wrad)
	  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, atypl, natyps,wrad)

	asa=asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc, atypstruc, natyps,wrad)
	grid.append(asar+asal-asa)
    np.save('surfacegrid-atoms-'+wout+data.split('.')[-2].strip('/')+'.npy',np.array(grid))
elif surface=='opls':
    mp = np.arange(0,83)
    replace={30:2, 31:3, 65:4, 66:5, 67:6, 68:7, 69:8, 70:9, 71:10, 80:11, 81:12, 82:13}
    mp[replace.keys()] = replace.values()
    modatypstruc = mp[modatypstruc]
    modatypl=mp[modatypl]
    modatypr=mp[modatypr]
    natyps=13
    grid=[]
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
      asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
	  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
	asa=asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc, modatypstruc, natyps,wrad)
	grid.append(asar+asal-asa)
    np.save('surfacegrid-'+surface+wout+'-'+data.split('.')[-2].strip('/')+'.npy',np.array(grid))
elif surface=='tobi':
    mp = np.arange(0,33)
    replace={32:19}
    mp[replace.keys()] = replace.values()
    modatypstruc = mp[modatypstruc]
    modatypl=mp[modatypl]
    modatypr=mp[modatypr]
    natyps=19
    grid=[]
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
      asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
	  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
	asa=asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc, modatypstruc, natyps,wrad)
	grid.append(asar+asal-asa)
    np.save('surfacegrid-'+surface+wout+'-'+data.split('.')[-2].strip('/')+'.npy',np.array(grid))
elif surface=='attract':
    natyps=32
    mp0 = np.arange(0,100)
    replace0={99:32}
    mp0[replace0.keys()] = replace0.values()
    modatypl = mp0[modatypl]
    modatypr = mp0[modatypr]
    modatypstruc = mp0[modatypstruc]
    cgname = ''
    if '--coarse_grained' in sys.argv:
      cgname = '-cg-radii'
      mp = np.arange(0.,33.)
      replace={1:2.0, 2:1.9, 3:1.95, 4:1.9, 5:1.9, 6:1.9, 7:2.0, 8:2.0, 9:1.9, 10:2.0, 11:1.9, 12:2.0, 13:1.9, 14:2.2, 15:2.2, 16:2.0, 17:1.9, 18:2.0, 19:2.0, 20:2.2, 21:2.2, 22:1.9, 23:1.9, 24:1.9, 25:2.0, 26:2.2, 27:2.2, 28:2.2, 29:2.0, 30:1.6, 31:1.5, 32:1.7}
      mp[replace.keys()] = replace.values()
      Rrvdw = mp[modatypr]
      Rlvdw = mp[modatypl]
      Rvdwstruc = mp[modatypstruc]

    grid=[]
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
      asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
	  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
	asa=asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc, modatypstruc, natyps,wrad)
	grid.append(asar+asal-asa)
    np.save('surfacegrid-'+surface+wout+'-'+data.split('.')[-2].strip('/')+cgname+'.npy',np.array(grid))
elif surface=='gaa':
    natyps=27
    grid=[]
    if modeon==False:
      coor = np.array(collectlib.collect_all_coor())
      coor_r=coor[:rlen]
      coor_l=coor[rlen:]
      asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
      asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
    while 1:
	if name is not None: 
          newargs = collectinsert + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
          if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
	    break
          collectlib.collect_iattract(newargs)
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	if modeon:
	  coor_r=coor[:rlen]
	  coor_l=coor[rlen:]
	  asar=asalib.asaatom(coor_r[:,0],coor_r[:,1],coor_r[:,2], Rrvdw, modatypr, natyps,wrad)
	  asal=asalib.asaatom(coor_l[:,0],coor_l[:,1],coor_l[:,2], Rlvdw, modatypl, natyps,wrad)
	asa=asalib.asaatom(coor[:,0],coor[:,1],coor[:,2], Rvdwstruc, modatypstruc, natyps,wrad)
	grid.append(asar+asal-asa)
    np.save('surfacegrid-'+surface+wout+'-'+data.split('.')[-2].strip('/')+'.npy',np.array(grid))
