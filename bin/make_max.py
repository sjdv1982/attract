#! /usr/bin/env python

import sys

h = sys.argv[1]
variables = {'MAXSTRUC':'','MAXRESTRAINTS':''}
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
	  if item in variables:
	    item = variables[item]
	  else:
            item = item.lower()
	  
      if item == '':
	output = ''
	break

      output += ' '+item
      
    if not output == '':
      print output
	  
	