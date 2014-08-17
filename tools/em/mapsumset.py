import sys, numpy
from em_iolib import read_situs, write_situs, write_mask

try:
  situsfile = sys.argv[1]
  targetfile = sys.argv[2]
  mapsum = float(sys.argv[3])
  data, gridspacing, origin = read_situs(situsfile)
except:
  import traceback
  traceback.print_exc()
  print
  print("Please provide SITUS density source map, target map, desired map sum")
  sys.exit()

curr_mapsum = numpy.sum(data, axis=None)
scale = mapsum / curr_mapsum
data *= scale
write_situs(targetfile, data, gridspacing, origin)