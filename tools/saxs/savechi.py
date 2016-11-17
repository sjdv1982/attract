import sys
import numpy

data = open(sys.argv[1]).readlines()
data = [float(l.split()[-1]) for l in data if 'Energy' in l]
data.sort()
data = numpy.array(data)
numpy.save(sys.argv[2],data)
if len(sys.argv) > 3:
  numpy.save(sys.argv[3],numpy.ones(len(data)))