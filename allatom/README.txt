This document describes the generation of ATTRACT parameters from CNS format.
See allatom.txt for more technical detail.


Input: 
- A CNS parameter file. This file must contain NONBONDED statements for all the desired atom types (other statements are ignored). 
Example: oplsx.par
- A CNS topology file. Example: topallhdg.pro

Steps:
1. Convert the parameter file to a trans file using the trans_par.py tool. The "start" parameter specifies an offset for ATTRACT atom types; 75 or 33 would make sense to avoid overlap with the standard reduced atom types.
Example: python trans_par.py oplsx.par 75
The trans file contains the ATTRACT atom type (a number), the VdW radius parameters, and a list of CNS atom types (code strings). Example: oplsx.trans
2. Edit the trans file. In general, it is best to assign different ATTRACT atom type numbers than 1-31 (the reduced ones). Example: oplsx.trans
3. Generate the ATTRACT parameters using the gen_par tool
Example: python gen_par.py oplsx.trans > allatom.par

4. Choose between 4a and 4b.
 4a. Instead of bin/reduce, use reduce.py to reduce your PDBs. This requires the trans file and the topology file.
Example: python reduce.py mypdb.pdb oplsx.trans topallhdg.pro 
This will create mypdb-aaX.pdb. 
The coordinates of all missing atoms are represented by XXXXXXXX. You will have to fill them in yourself.
 4b. Instead of bin/reduce, use pqreduce.py to reduce your PDBs. This requires the trans file and the topology file, and the installation of pdb2pqr.
Example: python pqreduce.py mypdb.pdb oplsx.trans topallhdg.pro 
This will create mypdb-aa.pdb. All missing will atoms have been filled in by pdb2pqr. 
5. Run ATTRACT normally, with the following differences:
- Use the full-atom "reduced" PDB files as generated above
- Instead of parmu.par, use allatom.par as generated above
- Use a constant dielectric by setting the --cdie flag. The recommended value is 10. This can be set with --epsilon 10.
- If you use a grid, use make-grid with --cdie and --epsilon 10 as well. It is also recommended to specify --alphabet <trans file> to save memory.
