import subprocess, os

subprocess.call('bash refin.sh',shell=True)
os.system('gawk -f diff.awk resultk.dat > result-amber.score')
for i in range(1000):
  subprocess.call('python rename_chain.py 2 bsf'+str(i+1)+'.pdb',shell=True)
  
os.system('python $ATTRACTTOOLS/joinmodel.py '+" ".join(['bsf'+str(i+1)+'-collect.pdb' for i in range(1000)])+' > result-amber.pdb')
import numpy as np
struc, score = np.loadtxt('result-amber.score',unpack=True,usecols=[0,1])
total = zip(score,struc)
total.sort()
os.system('python $ATTRACTTOOLS/joinmodel.py '+" ".join(['bsf'+str(int(i))+'-collect.pdb' for j,i in total])+' > result-amber-sorted.pdb')
out = open('result-amber-sorted.score')
for i, item in enumerate(total):
  out.write(str(i+1)+' '+str(item[0])+'\n')
  
out.close()
