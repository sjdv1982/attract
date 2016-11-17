#Performs IRMSD and fnat calculation on pepATTRACT AMBER refined docking models
pdb4amber -i receptor-refe.pdb -y -d -p -o murks.pdb
tleap -f leap.in
ambpdb -p exp.top < exp.crd > receptor-refe-amber.pdb
pdb4amber -i ligand-refe.pdb -y -d -p -o murks.pdb
tleap -f leap.in
ambpdb -p exp.top < exp.crd > ligand-refe-amber.pdb
python $ATTRACTTOOLS/irmsd-pdb.py result-amber-sorted.pdb receptor-refe-amber.pdb ligand-refe-amber.pdb > result-amber-sorted.irmsd
python $ATTRACTTOOLS/fnat-pdb.py result-amber-sorted.pdb 5 receptor-refe-amber.pdb ligand-refe-amber.pdb > result-amber-sorted.fnat