#! /usr/bin/env python

import sys

h = sys.argv[1]
variables = {'MAXSTRUC':'','MAXATOM':'maxatom','MAXRES':'maxres','TOTMAXATOM':'totmaxatom','TOTMAXRES':'totmaxres','MAXLIG':'maxlig','MAXMODE':'maxmode','MAXMOLPAIR':'maxmolpair','MAXDOF':'maxdof','MAXATOMTYPES':'maxatomtypes',
	     'MAXSELECTION':'maxselection','MAXRESTRAINTS':'','MAXENS':'maxens','MAXLENINDEXMODE':'maxlenindexmode','MAXINDEXMODE':'maxindexmode'}
for line in open(h).readlines():
  if 'typedef' in line:
    continue
  elif 'const int' in line:
    tmp = line.split()
    output = '      integer, parameter ::'
    for item in tmp:
      if '/' in item:
	break
      
      elif 'const' in item or 'int' in item:
	continue
      
      else:
	if ';' in item:
	  item = item.replace(';','')
	  
	if 'MAX' in item:
	  item = variables[item]
	  
      if item == '':
	output = ''
	break

      output += ' '+item
      
    if not output == '':
      print output
	  
	