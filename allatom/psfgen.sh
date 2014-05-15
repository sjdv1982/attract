#/bin/sh
grep 'ATOM' $1 >  inprotein.pdb
grep 'ATOM' $2 >  innucl.pdb
sed -i 's/OP1/O1P/g' innucl.pdb
sed -i 's/OP2/O2P/g' innucl.pdb	
vmd -dispdev text -e /home/ichauvot/attract/allatom/prot.pgn > vmdprotlog.out
vmd -dispdev text -e /home/ichauvot/attract/allatom/nucl.pgn > vmdnucllog.out
#rm inprotein.pdb
#rm innucl.pdb
