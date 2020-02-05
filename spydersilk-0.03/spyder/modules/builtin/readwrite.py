# Copyright 2011, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import os
from spyder.loader import file_load, file_access

def readCurrDir(filename, checkonly):
  okmods = os.F_OK + os.R_OK
  if file_access(filename, okmods) == False:
    return None
  if checkonly == True: return True
  try:
    return file_load(filename, "rb")
  except:
    return None

def readLocal(filename,checkonly):
  okmods = os.F_OK + os.R_OK
  if file_access(filename, okmods) == False:
    return None
  if checkonly == True: return True
  try:
    return file_load(filename, "rb")
  except:
    return None

if (spyder.python2):
  import urllib, urllib2
  urllib.request = spyder.newnamespace()
  urllib.request.urlopen = urllib2.urlopen
else:
  import urllib.request

def readRemote(filename, checkonly):
  try:
    f = urllib.request.urlopen(filename)
  except:
    return None
  if checkonly == True:
    f.close()
    return True
  return f
  
def writeCurrDir(filename, checkonly):
  e = os.path.exists(filename)
  try:
    f = open(filename, "ab")
    f.close()
    if checkonly == True:
      if not e: os.remove(filename)
      return True
    else: 
      f = open(filename, "wb")
      return f
  except:
    return None
  
def writeLocal(filename, checkonly):
  e = os.path.exists(filename)
  try:
    f = open(filename, "ab")
    f.close()
    if checkonly == True:
      if not e: os.remove(filename)
      return True
    else: 
      f = open(filename, "wb")
      return f
  except:
    return None
