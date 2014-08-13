import sys
from em_iolib import read_situs, write_situs
import numpy

packing = 0.824 #unit: Daltons per A**3


try:
  situsfile = sys.argv[1]
  targetfile = sys.argv[2]
  protein_mass = float(sys.argv[3])
  data, gridspacing, origin = read_situs(situsfile)
except:
  import traceback
  traceback.print_exc()
  print
  print("Please provide SITUS density source map, target map, protein mass")
  sys.exit()

voxelsize = gridspacing**3
negadens = numpy.amin(data)
print "SITUS voxel size: %.3f A**3" % voxelsize
if negadens < 0:
  print "Lowest voxel value: %.3f A**3" % negadens
proteinvolume = protein_mass / packing #Unit: A**3
print "Protein volume: %.3f A**3" % proteinvolume  
topvoxels =  int(proteinvolume/voxelsize+0.5)
print "Estimated protein voxels: %d" % topvoxels
data2 = numpy.sort(data, axis=None)
topsum = data2[-topvoxels:].sum()  
print "Sum of all voxels: %.3f" % data.sum()  
print "Sum of protein voxels: %.3f, rescaling to target protein mass..." % topsum
topsum += topvoxels * -negadens
factor = protein_mass/topsum
data = (data-negadens) * factor
write_situs(targetfile, data, gridspacing, origin)
