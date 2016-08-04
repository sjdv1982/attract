import os

os.system('python $ATTRACTDIR/../examples/attractsaxs/rescore.py result-run1 300.0')
os.system('python $ATTRACTDIR/../examples/attractsaxs/cluster.py result-run1-resorted-gaa4_6-joinedscore_300.0 6.5')