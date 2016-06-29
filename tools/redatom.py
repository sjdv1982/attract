#      program redatom
import sys     
#  run by:python redatom.py input.pdb >out.pdb

if '--atomtypefile' in sys.argv:
  atomtypes = sys.argv[sys.argv.index('--atomtypefile')+1]
else:
  print 'please give file with all the atomtypes for the model'
  sys.exit()
  
build_atom = open(atomtypes, 'r')
f1lines=build_atom.readlines()

if '--output' in sys.argv:
  outname=sys.argv[sys.argv.index('--output')+1]
else:
  print 'give outputname'
  sys.exit()

if '--pdb' in sys.argv:
  eingabe=sys.argv[sys.argv.index('--pdb')+1]
else: 
  print 'please give pdb to convert'
  sys.exit()
  
  
  
if '--addNH' in sys.argv:
  f=open(eingabe, 'r') 
  i=0
  iflater=0
  flines=[]
  for line in f:
      if (line[:4]=='ATOM') and line[13:16]!='OXT' and line[13:14]!='H' and line[12:13]!='H':
	  flines.append(line)
	  i+=1

  f.close()
  # automatic addition of hydrogens
  # treat first N separately
  natom=i
  nnew=0
  for a in range(natom):
      if flines[a][13:16]=='N  ' and flines[a][17:20]!='PRO':
	  nnew+=1
  for a in range(natom+nnew):
      if flines[a][13:16]=='N  ' and flines[a][17:20]!='PRO':
	  flines.insert(a+1,flines[a][:13]+'%-3s' % 'H'+flines[a][16:30]+'%8s' % round(float(flines[a][30:38])*
			1.66-float(flines[a+1][30:38])*0.66,3)+'%8s' % round(float(flines[a][38:46])*1.66-float(flines[a+1][38:46])*0.66,3)+
			'%8s' % round(float(flines[a][46:54])*1.66-float(flines[a+1][46:54])*0.66,3)+flines[a][54:69]+'\n')

else:
  fobj=open(eingabe, 'r')
  lines=fobj.readlines()
  i=0
  flines=[]
  for line in lines:
    if line[:4]=='ATOM':
      flines.append(line)
      i+=1
  nnew=0
  natom=i
  fobj.close()
#
#     generierung der pdbcade tabelle und der integercode tabellen
#     fuer jede Aminosaeure;
#
if '--atomtypesonly' in sys.argv:
  outobj=open(outname,'w+')
  t=0
  atm=natom
  test=0
  for i in range(natom):
      for j in range(len(f1lines)):
	  if f1lines[j][:4]==flines[i][17:21] and f1lines[j][6:10]==flines[i][12:16]:
	      test+=1
	      outobj.write(flines[i][:55]+f1lines[j][12:16]+flines[i][59:67]+" "+"0"+" "+"1.00"+"\n")

      if flines[i][:3]=='TER' or flines[i][:3]=='END':
	  outobj.write(flines[i])
	  t=0
	  atm=0
	  test=0

  outobj.close()
else:  
  outobj=open(outname,'w+')
  t=0
  atm=0
  test=0
  for i in range(natom+nnew):
      atm+=1

      if i<(natom+nnew-1):
	  if flines[i][12:16]==' N  ':
	      t+=1
      for j in range(len(f1lines)):
	  if f1lines[j][:4]==flines[i][17:21] and f1lines[j][6:10]==flines[i][12:16]:
	      test+=1
	      outobj.write(flines[i][:6]+'%5s' % atm+flines[i][11:21]+'%5s' % t+flines[i][26:55]+f1lines[j][12:16]+f1lines[j][17:25]+" "+"0"+" "+"1.00"+"\n")

      if flines[i][:3]=='TER' or flines[i][:3]=='END':
	  outobj.write(flines[i])
	  t=0
	  atm=0
	  test=0

  outobj.close()
if test!=atm:
  print 'something went wrong, wrote: ', test, atm
	