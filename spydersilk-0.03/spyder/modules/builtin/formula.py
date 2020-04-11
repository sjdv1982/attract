# Copyright 2008-2011, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

from math import *

def formula(s, params):
  p = ",".join(params)
  spyder.safe_eval.safe_eval(p)
  spyder.safe_eval.safe_eval(s)  
  funcdef = """def func(%s):
  return %s
""" % (p, s)
  exec(funcdef)
  return func
  

