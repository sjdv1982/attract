def generate_narefine(m):
  import os
  pattern = os.path.splitext(m.pdbfile.name)[0]
  ret = """#!/bin/bash -i
set -u -e

fixprotein=$1
d='.'
pattern=%s
nparallel=%s
use_cuda=%s
use_5PO=%s
""" % (pattern, m.nparallel, int(m.use_cuda), int(m.use_5PO))
  ret += """
nmodels=`awk '$1 == "MODEL"' $d/$pattern.pdb | wc -l`

echo '**************************************************************'
echo 'Prepare files...'
echo '**************************************************************'
if use_5PO=0; then
  gawk -f $NAREFINE/scripts/withdraw5PO.awk $d/$pattern.pdb > $d/$pattern-no5PO.pdb
  pattern=$pattern-no5PO  
fi  
$NAREFINE/prep.sh $d $d/$pattern.pdb $use_5PO $d/$pattern.top $d/$pattern $nparallel $fixprotein
python $NAREFINE/scripts/generate_infiles.py $d/mini1.in $d/mini2.in %.1f %d %d $d/md1-X.in $d/md2-X.in %.1f %.1f %s %s %d %d %.1f %.1f %d $d/minc.in $d/mingb.in  
""" % (m.offset, m.min2_cycles, int(m.use_simulation),
       m.md1.restraints_value, m.md2.restraints_value, 
       m.md1.restraints_atoms, m.md2.restraints_atoms,
       m.md1.nsteps, m.md2.nsteps,
       m.md1.temperature, m.md2.temperature, 
       m.minc_cycles
      )
  ret += """
cd $d

echo '**************************************************************'
echo 'First minimization...'
echo '**************************************************************'
$NAREFINE/mini1.sh $nmodels $pattern.top $pattern $nparallel

echo '**************************************************************'
echo 'Second minimization...'
echo '**************************************************************'
$NAREFINE/mini2.sh $nmodels $pattern.top $pattern $nparallel $use_cuda

echo '**************************************************************'
echo 'Scoring of the minimization...'
echo '**************************************************************'
$NAREFINE/score-mini.sh $nmodels $pattern.top $pattern $nparallel $use_5PO  

echo '**************************************************************'
echo 'MD simulation...'
echo '**************************************************************'
$NAREFINE/md.sh $nmodels $pattern.top $pattern $nparallel $use_cuda

echo '**************************************************************'
echo 'Scoring of the final results...'
echo '**************************************************************'
$NAREFINE/score.sh $nmodels $pattern.top $pattern $nparallel $use_5PO
  """
  return ret