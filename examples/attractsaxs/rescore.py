import os, subprocess, sys
import numpy as np
name = sys.argv[1]
weight = float(sys.argv[2])
workdir = os.getcwd()
os.environ["PYTHONPATH"]=os.environ["PYTHONPATH"]+':'+os.environ["ATTRACTDIR"]+':'+os.environ["ATTRACTTOOLS"]+':'+os.environ["ATTRACTDIR"]+'/../examples/attractsaxs/LX_Scoring/Rescoring'
if not os.path.exists('ubA-gaa.pdb'):
  os.system('python $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/Grids/redatom.py --pdb ubB-aa.pdb --atomtypefile $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/Grids/atomtypes_ga.dat --output ubB-gaa.pdb')
  os.system('python $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/Grids/redatom.py --pdb ubA-aa.pdb --atomtypefile $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/Grids/atomtypes_ga.dat --output ubA-gaa.pdb')
  
subprocess.call('python $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/Rescoring/@rank.py --input '+name+'.dat --proteinmodel gaa ubA-gaa.pdb ubB-gaa.pdb --steppotential $ATTRACTDIR/../examples/attractsaxs/LX_Scoring/params/MC_gaa_6_10.par --bins 2 4 6  --rescore --output '+name+'-gaa4_6.rescore',shell=True)
score1 = open(name+'.dat').readlines()
score1 = [float(l.split()[-1]) for l in score1 if 'Energy' in l]
score1 = np.array(score1)
score2 = open(name+'-gaa4_6.rescore').readlines()
score2 = [float(l.split()[-1]) for l in score2 if 'Energy' in l]
score2 = np.array(score2)
ene = weight*score1+score2

out = open(name+'-resorted-gaa4_6-joinedscore_'+str(weight)+'.rescore','w')
for e in ene:
  out.write('Energy: '+str(e)+'\n')
  
out.close()
os.system('python $ATTRACTTOOLS/fill-energies.py '+name+'.dat '+name+'-resorted-gaa4_6-joinedscore_'+str(weight)+'.rescore > '+name+'-resorted-gaa4_6-rescored-joinedscore_'+str(weight)+'.dat')
os.system('python $ATTRACTTOOLS/sort.py '+name+'-resorted-gaa4_6-rescored-joinedscore_'+str(weight)+'.dat > '+name+'-resorted-gaa4_6-joinedscore_'+str(weight)+'.dat')