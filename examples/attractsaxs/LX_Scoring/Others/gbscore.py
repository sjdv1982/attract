#collect.py
import numpy as np
import numpy.ma as ma
import sys
from scipy.spatial.distance import cdist
import collectlibpy as collectlib
import gbenbeta

receptor=sys.argv[1]
ligand=sys.argv[2]
data=sys.argv[3]

collectlib.collect_init([data,receptor,ligand])

Sr=[]
Rrvdw=[]
charger=[]
#test=[]
fobj=open(receptor, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
#    test.append([float(fline[30:38]),float(fline[38:46]),float(fline[46:54])])
    charger.append(float(fline[60:67]))
    if fline[13]=='N':
      Rrvdw.append(1.6)
      Sr.append(0.76)
    elif fline[13]=='C':
      Rrvdw.append(1.7)
      Sr.append(0.76)
    elif fline[13]=='O':
      Rrvdw.append(1.5)
      Sr.append(0.83)
    elif fline[13]=='S':
      Rrvdw.append(1.9)
      Sr.append(0.96)
    elif fline[13]=='H' or fline[12]=='H':
      Rrvdw.append(1.2)
      Sr.append(0.8)
    elif fline[13]=='P':
      Rrvdw.append(2.0)
      Sr.append(0.86)
    else:
      print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', receptor
      sys.exit()
fobj.close()
#test=np.array(test)
Sr=np.array(Sr)
Rrvdw=np.array(Rrvdw)
charger=np.array(charger)

Sl=[]
Rlvdw=[]
chargel=[]
fobj=open(ligand, 'r')
for i,fline in enumerate(fobj):
  if fline[:4]=='ATOM':
    chargel.append(float(fline[60:67]))
    if fline[13]=='N':
      Rlvdw.append(1.6)
      Sl.append(0.76)
    elif fline[13]=='C':
      Rlvdw.append(1.7)
      Sl.append(0.76)
    elif fline[13]=='O':
      Rlvdw.append(1.5)
      Sl.append(0.83)
    elif fline[13]=='S':
      Rlvdw.append(1.9)
      Sl.append(0.96)
    elif fline[13]=='H' or fline[12]=='H':
      Rlvdw.append(1.2)
      Sl.append(0.8)
    elif fline[13]=='P':
      Rlvdw.append(2.0)
      Sl.append(0.86)
    else:
      print fline[13], ' is a new detected atomtype!!! in line ', i, 'of ', ligand
      sys.exit()
fobj.close()
Sl=np.array(Sl)
Rlvdw=np.array(Rlvdw)
chargel=np.array(chargel)

if '--Smulti' in sys.argv:
  multi = float(sys.argv[sys.argv.index('--Smulti')+1])
  Sl = Sl*multi
  Sr = Sr*multi
elif '--Snone' in sys.argv:
  multi = float(sys.argv[sys.argv.index('--Snone')+1])
  Sl = np.ones(np.shape(Sl))*multi
  Sr = np.ones(np.shape(Sr))*multi
  
  
Rvdwstruc=np.append(Rrvdw,Rlvdw)
Sstruc=np.append(Sr,Sl)
chargestruc=np.append(charger,chargel)
rlen = len(charger)
llen = len(chargel)
nstruc=0

if '--delta' in sys.argv:
  delta = float(sys.argv[sys.argv.index('--delta')+1])
else: 
  delta = 0.12
  

if '--gbscore' in sys.argv:
  if '--brfixed' in sys.argv:
    br = float(sys.argv[sys.argv.index('--brfixed')+1])
#    coor = np.array(collectlib.collect_all_coor())
#    coor_r=coor[:rlen]
#    coor_l=coor[rlen:]
#    enr=gbenbeta.gbenergy(coor_r, np.ones(rlen)*2.5, charger)
#    enl=gbenbeta.gbenergy(coor_l, np.ones(llen)*2.5, chargel)
    while 1:
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	coor_r = coor[:rlen]
	coor_l = coor[rlen:]
	en=gbenbeta.gbenfixed(coor_r, coor_l, br, charger, chargel)
	#en2=gbenbeta.gbenergy(coor, np.ones(rlen+llen)*2.5, chargestruc)-enr-enl
	print en #, en2
  else:
    coor = np.array(collectlib.collect_all_coor())
    coor_r=coor[:rlen]
    coor_l=coor[rlen:]
    bradr=gbenbeta.gbrad(coor_r, Rrvdw, Sr, delta)
    #for i in range(len(bradr)):
    #  print i,bradr[i]
    #sys.exit()
    bradl=gbenbeta.gbrad(coor_l, Rlvdw, Sl, delta)
    enr=gbenbeta.gbenergy(coor_r, bradr, charger)
    enl=gbenbeta.gbenergy(coor_l, bradl, chargel)
    #print 'gbreceptor ', enr
    #print 'gbligand ',enl
    #bradstruc=np.append(bradr,bradl)
    while 1:
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	bradstruc=gbenbeta.gbradcomplex(coor, Rrvdw, Rlvdw, Sr, Sl, bradr, bradl, delta) #exact calculation
	en=gbenbeta.gbenergy(coor, bradstruc, chargestruc)
	#en=gbenbeta.gbencmpx(coor, rlen, Rvdwstruc, Sstruc, chargestruc, bradstruc)
	print enl+enr-en

if '--electrostatic' in sys.argv:
    if '--shift' in sys.argv:
      shift=float(sys.argv[sys.argv.index('--shift')+1])
    else:
      shift=0.
    #coor = np.array(collectlib.collect_all_coor())
    #coor_r=coor[:rlen]
    #coor_l=coor[rlen:]
    #enr=gbenbeta.gbenergy(coor_r, np.zeros(rlen), charger)
    #enl=gbenbeta.gbenergy(coor_l, np.zeros(llen), chargel)
    while 1:
	result = collectlib.collect_next()
	if result: break
	nstruc+=1
	coor = np.array(collectlib.collect_all_coor())
	coor_r = coor[:rlen]
	coor_l = coor[rlen:]
	en=gbenbeta.elec(coor_r, coor_l, charger, chargel, shift)
	#en2=gbenbeta.gbenergy(coor, np.zeros(rlen+llen), chargestruc)-enr-enl
	print en
	
	
	
