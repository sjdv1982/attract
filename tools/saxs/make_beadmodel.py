import os
import sys
from multiprocessing import Pool

def run_dammif(indata):
  i, args = indata
  os.system(atsasdir+'/dammif --prefix=model-'+str(i+1)+' --mode=fast pr.out '+args)
  
saxsdata = sys.argv[1]
args = ''
if len (sys.argv)> 2:
  args = " ".join(sys.argv[2:])

atsasdir = os.environ["ATSASDIR"]
if atsasdir == '' or not os.path.exists(atsasdir+'/datgnom'):
  print "ERROR: Please ensure that you have the ATSAS software installed and the environment variable ATSASDIR set correctly!"
  sys.exit(1)

os.system('mkdir dammifrun')
os.chdir('dammifrun')
os.system(atsasdir+'/datgnom4 '+saxsdata+' -o pr.out')
#run dammif 20 times

p = Pool(8)
indata = [ (i,args) for i in range(20)]
p.map(run_dammif,indata)
p.close()
p.join()
os.system(atsasdir+'/damaver model*-1.pdb --automatic '+args)
os.chdir('..')
os.system('ln -s dammifrun/damaver.pdb beadmodel.pdb')

