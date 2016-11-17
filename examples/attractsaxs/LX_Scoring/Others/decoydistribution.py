#decoydistribution.py
import sys,os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

x=sys.argv[1]
xvalue=sys.argv[2]
border=int(sys.argv[3])
structure=sys.argv[4]

dirlist = [z for z in os.listdir('.') if os.path.isdir(z)]
#dirlist.remove('1BGX')
#dirlist.remove('2HMI')
dirlist = np.array(dirlist, dtype=str)

classon = False
if '--classification' in sys.argv:
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
    classes = np.zeros(len(dirlist),dtype = int)
    dictio = { "rigid":0, "medium":1, "hard":2}
    for i, folder in enumerate(dirlist):
      classes[i] = dictio[comclass[1][comclass[0].index(folder)]]
  elif classtype == 'proteintype':
    setname=['enzyme', 'antibody', 'other']
    classes = np.zeros(len(dirlist),dtype=int)
    dictio = { "enzyme":0, "antibody":1, "other":2}
    for i, folder in enumerate(dirlist):
      classes[i] = dictio[comclass[2][comclass[0].index(folder)]]
  dirlist = np.array(dirlist)
  sort = np.argsort(classes)
  dirlist = dirlist[sort]
  classlabel = np.array(setname,dtype = str)
  classlabel = classlabel[classes[sort]]
  trainlen = len(ma.compressed(ma.masked_not_equal(classes,0)))
  testlen = len(ma.compressed(ma.masked_not_equal(classes,1)))
  otherlen = len(ma.compressed(ma.masked_not_equal(classes,2)))

print setname[0],trainlen
print setname[1],testlen
print setname[2],otherlen

if '--plotdensity' in sys.argv:
  if structure!='benchmark':
      if xvalue=='irmsd' or xvalue=='lrmsd':
	xi=np.genfromtxt(x)[:border,-1]
      else:
	xi=np.genfromtxt(x)[:border]
	
      fig=plt.figure()
      ax=fig.add_subplot(111)
      n, bins, patches = ax.hist(xi, 100, normed=False, facecolor='green', alpha=0.75)
      ax.set_xlabel(xvalue)
      ax.set_ylabel('density')

      plt.show()

  else:
      xi=[]
      yi=[]
      good=[]
      for direct in dirlist:
	  if xvalue=='irmsd' or xvalue=='lrmsd':
	    xii=np.zeros(102)
	    yii=np.zeros(102)
	    yii[1:-1],xii[1:]=np.histogram(np.genfromtxt(direct+'/'+x)[:border,-1],bins=100)
	    xi.append(xii)
	    yi.append(yii)
	  else:
	    xii=np.zeros(102)
	    yii=np.zeros(102)
	    yii[1:-1],xii[1:]=np.histogram(np.genfromtxt(direct+'/'+x)[:border],bins=100)
	    xi.append(xii)
	    yi.append(yii)
      xi=np.array(xi)
      yi=np.array(yi)
      

      
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.autoscale(True)
      frame = 0
      ax.set_title(dirlist[frame])
      ax.grid(True)
      ax.set_xlabel(xvalue)
      #ax.set_ylim([-100,20])
      plt.subplots_adjust(left=0.25, bottom=0.25)
      ln,=ax.step(xi[frame],yi[frame],where='post')
      axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
      sframe = Slider(axframe, 'Frame', 0, len(dirlist)-1, valinit=0,valfmt='%d')
	      # call back function
      def update(val):
	  frame = np.floor(sframe.val)
	  ln.set_xdata(xi[frame])
	  ln.set_ydata(yi[frame])
	  ax.set_title(dirlist[frame])
	  ax.grid(True)
	  ax.autoscale(True)
	  ax.relim()
	  ax.autoscale_view()
	  #ax.set_ylim([-100,20])
	  plt.draw()
	      # connect callback to slider   
      sframe.on_changed(update)
      plt.show()

   
if '--statistics' in sys.argv:
    if structure=='benchmark':
      if xvalue=='irmsd':
	if classon:
	  print xvalue,'|  class |   < 1|   < 2|   < 3|   < 4|   < 5|   > 5| total|'
	else:
	  print xvalue,'|   < 1|   < 2|   < 3|   < 4|   < 5|   > 5| total|'
	
	for i, direct in enumerate(dirlist):
	  xi=np.genfromtxt(direct+'/'+x)[:border,-1]
	  num1=len(ma.compressed(ma.masked_outside(xi,0,1)))
	  num2=len(ma.compressed(ma.masked_outside(xi,1,2)))
	  num3=len(ma.compressed(ma.masked_outside(xi,2,3)))
	  num4=len(ma.compressed(ma.masked_outside(xi,3,4)))
	  num5=len(ma.compressed(ma.masked_outside(xi,4,5)))
	  total=len(xi)
	  test=len(ma.compressed(ma.masked_less_equal(xi,5)))
	  if classon:
	    print direct, "%10s" % classlabel[i], "%6s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%6s" % str(num4),"%6s" % str(num5), "%6s" % str(test), "%6s" % str(total)
	  else:
	    print direct, "%8s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%6s" % str(num4),"%6s" % str(num5), "%6s" % str(test), "%6s" % str(total)
      elif xvalue=='lrmsd':
	if classon:
	  print xvalue,'|  class |   < 2|   < 4|   < 6|   < 8|   < 10|   > 10| total|'
	else:  
	  print xvalue,'|   < 2|   < 4|   < 6|   < 8|   < 10|   > 10| total|'
	for i, direct in enumerate(dirlist):
	  xi=np.genfromtxt(direct+'/'+x)[:border,-1]
	  num1=len(ma.compressed(ma.masked_outside(xi,0,2)))
	  num2=len(ma.compressed(ma.masked_outside(xi,2,4)))
	  num3=len(ma.compressed(ma.masked_outside(xi,4,6)))
	  num4=len(ma.compressed(ma.masked_outside(xi,6,8)))
	  num5=len(ma.compressed(ma.masked_outside(xi,8,10)))
	  total=len(xi)
	  test=len(ma.compressed(ma.masked_less_equal(xi,10)))
	  if classon:
	    print direct, "%10s" % classlabel[i], "%6s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%7s" % str(num4),"%7s" % str(num5), "%6s" % str(test), "%6s" % str(total)
	  else:
	    print direct, "%8s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%7s" % str(num4),"%7s" % str(num5), "%6s" % str(test), "%6s" % str(total)

      elif xvalue=='capristars':
	if classon:
	  print xvalue,'|  class | high | medium | acceptable | notacceptable | total|'
	else:
	  print xvalue,'| high | medium | acceptable | notacceptable | total|'
	for i,direct in enumerate(dirlist):
	  xi=np.genfromtxt(direct+'/'+x)[:border]
	  num1=len(ma.compressed(ma.masked_not_equal(xi,3)))
	  num2=len(ma.compressed(ma.masked_not_equal(xi,2)))
	  num3=len(ma.compressed(ma.masked_not_equal(xi,1)))
	  num4=len(ma.compressed(ma.masked_not_equal(xi,0)))
	  total=len(xi)
	  if classon:
	    print direct, "%10s" % classlabel[i], "%10s" % str(num1),"%8s" % str(num2),"%12s" % str(num3),"%15s" % str(num4),"%7s" % str(total)
	  else:
	    print direct, "%12s" % str(num1),"%8s" % str(num2),"%12s" % str(num3),"%15s" % str(num4),"%7s" % str(total)
	    
      elif xvalue=='fnat':
	if classon:
	  print xvalue,'|  class |  > 0.8|  > 0.6|  > 0.4|  > 0.3|  > 0.2|  > 0.1|  < 0.1| total|'
	else:
	  print xvalue,'|  > 0.8|  > 0.6|  > 0.4|  > 0.3|  > 0.2|  > 0.1|  < 0.1| total|'
	  
	for i, direct in enumerate(dirlist):
	  xi=np.genfromtxt(direct+'/'+x)[:border]
	  num1=len(ma.compressed(ma.masked_outside(xi,0.1,0.2000001)))
	  num2=len(ma.compressed(ma.masked_outside(xi,0.2,0.3000001)))
	  num2b=len(ma.compressed(ma.masked_outside(xi,0.3,0.4000001)))
	  num3=len(ma.compressed(ma.masked_outside(xi,0.4,0.60000001)))
	  num4=len(ma.compressed(ma.masked_outside(xi,0.6,0.80000001)))
	  num5=len(ma.compressed(ma.masked_outside(xi,0.8,1.00000001)))
	  total=len(xi)
	  test=len(ma.compressed(ma.masked_greater(xi,0.1)))
	  if classon:
	    print direct, "%10s" % classlabel[i], "%6s" % str(num5),"%7s" % str(num4),"%7s" % str(num3),"%7s" % str(num2b),"%7s" % str(num2),"%7s" % str(num1), "%8s" % str(test), "%6s" % str(total)
	  else:
	    print direct, "%8s" % str(num5),"%7s" % str(num4),"%7s" % str(num3),"%7s" % str(num2b),"%7s" % str(num2),"%7s" % str(num1), "%8s" % str(test), "%6s" % str(total)



    else:
      if xvalue=='irmsd':
	xi=np.genfromtxt(x)[:border,-1]
	num1=len(ma.compressed(ma.masked_outside(xi,0,1)))
	num2=len(ma.compressed(ma.masked_outside(xi,1,2)))
	num3=len(ma.compressed(ma.masked_outside(xi,2,3)))
	num4=len(ma.compressed(ma.masked_outside(xi,3,4)))
	num5=len(ma.compressed(ma.masked_outside(xi,4,5)))
	total=len(xi)
	test=len(ma.compressed(ma.masked_less_equal(xi,5)))
	print xvalue,'|   < 1|   < 2|   < 3|   < 4|   < 5|   > 5| total|'
	print structure, "%8s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%6s" % str(num4),"%6s" % str(num5), "%6s" % str(test), "%6s" % str(total)
      elif xvalue=='lrmsd':
	xi=np.genfromtxt(x)[:border,-1]
	num1=len(ma.compressed(ma.masked_outside(xi,0,2)))
	num2=len(ma.compressed(ma.masked_outside(xi,2,4)))
	num3=len(ma.compressed(ma.masked_outside(xi,4,6)))
	num4=len(ma.compressed(ma.masked_outside(xi,6,8)))
	num5=len(ma.compressed(ma.masked_outside(xi,8,10)))
	total=len(xi)
	test=len(ma.compressed(ma.masked_less_equal(xi,10)))
	print xvalue,'|   < 2|   < 4|   < 6|   < 8|   < 10|   > 10| total|'
	print structure, "%8s" % str(num1),"%6s" % str(num2),"%6s" % str(num3),"%7s" % str(num4),"%7s" % str(num5), "%6s" % str(test), "%6s" % str(total)
      elif xvalue=='capristars':
	xi=np.genfromtxt(x)[:border]
	num1=len(ma.compressed(ma.masked_not_equal(xi,3)))
	num2=len(ma.compressed(ma.masked_not_equal(xi,2)))
	num3=len(ma.compressed(ma.masked_not_equal(xi,1)))
	num4=len(ma.compressed(ma.masked_not_equal(xi,0)))
	total=len(xi)
	print xvalue,'| high | medium | acceptable | notacceptable | total|'
	print structure, "%12s" % str(num1),"%8s" % str(num2),"%12s" % str(num3),"%15s" % str(num4),"%7s" % str(total)
      elif xvalue=='fnat':
	xi=np.genfromtxt(x)[:border]
	num1=len(ma.compressed(ma.masked_outside(xi,0.1,0.2)))
	num2=len(ma.compressed(ma.masked_outside(xi,0.2,0.4)))
	num3=len(ma.compressed(ma.masked_outside(xi,0.4,0.6)))
	num4=len(ma.compressed(ma.masked_outside(xi,0.6,0.8)))
	num5=len(ma.compressed(ma.masked_outside(xi,0.8,1.0)))
	total=len(xi)
	test=len(ma.compressed(ma.masked_greater(xi,0.1)))
	print xvalue,'|   > 0.8|   > 0.6|   > 0.4|   > 0.2|   > 0.1|   < 0.1| total|'
	print structure, "%9s" % str(num5),"%8s" % str(num4),"%8s" % str(num3),"%8s" % str(num2),"%8s" % str(num1), "%8s" % str(test), "%6s" % str(total)



















