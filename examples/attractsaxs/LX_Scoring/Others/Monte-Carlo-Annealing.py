import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy as np
import sys
from matplotlib import cm

datfile=sys.argv[1]

scores=np.genfromtxt(datfile, skip_header=1)
obj=open(datfile, 'r')
names=obj.readline()
names=names.strip()
names=names.split()

obj.close()

length=len(scores[0])
trainlength=(length/5)*4
testlength=length-trainlength

average=np.mean(scores[:,:trainlength],axis=1)
crossav=np.mean(scores[:,trainlength:],axis=1)

fig=plt.figure('Trainingset')
ax=fig.add_subplot(111)
ax.set_xlabel('Step')
ax.set_ylabel('testscore')
p=[ax.plot(scores[:,i],color=cm.hsv(1.*(i)/float(trainlength)), linewidth=1.5, label=names[i]) for i in range(trainlength)]
ax.legend(loc="upper left", bbox_to_anchor=(1,1), ncol=3, prop={'size':10})
ax.grid()

fig1=plt.figure('Trainingset-average')
ax1=fig1.add_subplot(111)
ax1.set_xlabel('Step')
ax1.set_ylabel('average testscore')
p1=ax1.plot(average)
ax1.grid()

fig2=plt.figure('Crossvalidation')
ax2=fig2.add_subplot(111)
ax2.set_xlabel('Step')
ax2.set_ylabel('testscore')
p2=[ax2.plot(scores[:,i],color=cm.hsv(1.*(i-trainlength)/float(testlength/2)), linewidth=1.5, label=names[i]) for i in range(trainlength,trainlength+testlength/2)]
ax2.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
ax2.grid()

fig4=plt.figure('Crossvalidation2')
ax4=fig4.add_subplot(111)
ax4.set_xlabel('Step')
ax4.set_ylabel('testscore')
p4=[ax4.plot(scores[:,i],color=cm.hsv(1.*(i-trainlength-testlength/2)/float(testlength/2)), linewidth=1.5, label=names[i]) for i in range(trainlength+testlength/2,length)]
ax4.legend(loc="upper left", bbox_to_anchor=(1,1), prop={'size':10})
ax4.grid()

fig3=plt.figure('Crossvalidation-average')
ax3=fig3.add_subplot(111)
ax3.set_xlabel('Step')
ax3.set_ylabel('average testscore')
p3=ax3.plot(crossav)
ax3.grid()

plt.show()