This document describes the generation of ATTRACT topology and parameters from CNS or OpenMM format.
See allatom.txt for more technical detail.

*********************************************************************************************************
Input:
CNS:
	- A CNS parameter file. This file must contain NONBONDED statements for all the desired atom types (other statements are ignored).
	Example: oplsx.par
	- A CNS topology file. Example: topallhdg.pro
OR: OpenMM:
	- An OpenMM .ffxml file. This file must contain a <NonbondedForce> section with <Atom> tags.
	  All <Atom> tags must contain a "class" or a "type" field. If there are "type" fields, the .ffxml file must also contain
		 an <AtomTypes> section that maps the type to a class.
		All atom classes must have the same Lennard-Jones parameters.

Steps:

1. Convert the parameter file to a trans file using the trans_par.py (CNS) or trans_par-ffxml (OpenMM TODO) tool.
	The "start" parameter specifies an offset for ATTRACT atom types; 75 or 33 would make sense to avoid overlap with the standard reduced atom types.
	Example: python trans_par.py oplsx.par 75
	The trans file contains the ATTRACT atom type (a number), the VdW radius parameters, and a list of CNS atom types (code strings). Example: oplsx.trans

2. Edit the trans file. In general, it is best to assign different ATTRACT atom type numbers than 1-31 (the reduced ones). Example: oplsx.trans
To create a combined file for DNA, RNA and protein, concatenate different trans files. Example: allatom.trans (from opls.trans and dna-rna.trans)

3. Generate the ATTRACT parameters using the gen_par tool
	Example: python gen_par.py allatom.trans > allatom.par



C. Usage
=====================================================================================================================================

4. Use allatom/aareduce.py to reduce your PDBs. This requires the trans file and the ATTRACT topology file.
	 python aareduce.py mypdb.pdb --trans oplsx.trans --top topallhdg5.3_HIS.pro
	 This will create mypdb-aaX.pdb.

   For reducing rna (resp. dna), use the topology files rna.top (resp. dna.top).
   The difference between those files consist in the charge of the sugar atom C2' (type C2R (resp. C2D)).

	 OpenMM example:

	 TODO: Not yet implemented in aareduce.py

	 Missing atoms:
	 Missing atoms give an error, unless --manual is specified.
   In that case, the coordinates of all missing atoms are represented by XXXXXXXX. You will have to fill them in yourself.
   Alternatively, you can use a backend (pdb2pqr or the WHATIF server) to complete hydrogens and other missing atoms


5. Run ATTRACT normally, with the following differences:
	- Use the all-atom "reduced" PDB files as generated above
	- Instead of parmw.par, use allatom.par as generated above, or dna-rna-prot-allatom.par
	- Use a constant dielectric by setting the --cdie flag. The recommended value is 10. This can be set with --epsilon 10.
	- If you use a grid, use make-grid with --cdie and --epsilon 10 as well. It is also recommended to specify --alphabet <trans file> to save time and space.

*********************************************************************************************************
