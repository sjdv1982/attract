python $ATTRACTTOOLS/saxs/make_weight.py exp-ori.dat > exp-processed.dat
python $ATTRACTTOOLS/saxs/make_beadmodel.py ../exp-processed.dat
python $ATTRACTDIR/../allatom/aareduce.py ubA.pdb --chain A --pdb2pqr --dumppatch > ubA.mapping
python $ATTRACTDIR/../allatom/aareduce.py ubB.pdb --chain B --pdb2pqr --dumppatch > ubB.mapping