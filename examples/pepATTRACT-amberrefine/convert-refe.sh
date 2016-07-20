pdb4amber -i receptor-refe.pdb -y -d -p -o murks.pdb
tleap -f leap.in
ambpdb -p exp.top < exp.crd > receptor-refe-amber.pdb
pdb4amber -i ligand-refe.pdb -y -d -p -o murks.pdb
tleap -f leap.in
ambpdb -p exp.top < exp.crd > ligand-refe-amber.pdb
