import sys,os
import random
import math

args = ' '.join(sys.argv[3:-1])
outfile = sys.argv[-1]

k = random.randint(999,99999)
jobsize = 25000
data = open(sys.argv[2]).readlines()
data =[l for l in data if 'Energy' in l]
chunks = int(math.ceil(len(data)/float(jobsize)))
os.system('python2 $ATTRACTTOOLS/split.py '+sys.argv[2]+' tmp'+str(k)+' '+str(chunks))
for i in range(1,chunks+1):
  os.system('python2 $ATTRACTTOOLS/saxs/saxs-score.py '+sys.argv[1]+' tmp'+str(k)+'-'+str(i)+' '+args+' | awk \'{print "Energy:", $2}\' > tmpene'+str(k)+'-'+str(i))

os.system('python2 $ATTRACTTOOLS/join.py tmpene'+str(k)+' --score > '+outfile)
os.system('rm tmp*'+str(k)+'*')
