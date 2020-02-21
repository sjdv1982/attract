# Copyright 2008, 2009, 2011 Sjoerd de Vries
# This file is part of the Spyder module: "atom" 
# For licensing information, see LICENSE.txt 

"""
Python script to contact the WHATIF web servers at the CMBI
In principle, this will work with any WHATIF server that accepts a single file and outputs a single file

Sjoerd de Vries, December 2006, ADAPTED November 2007 for the HADDOCK server

syntax: whatif-server.py <PDB-file> <WHATIF request> <WHATIF output name>

List of WHATIF server requests and output names:
    Add missing hydrogens: 
      request: htopo
      output name: hadded.pdb
    Completing sidechains:
      request: corall
      output name: fixed.pdb
      
    
To find out the WHATIF request for the server that you need: 
  - Go to http://swift.cmbi.ru.nl/WIWWWI/ and click on your server 
  - Click on View Frame Source
  - Look for this line:
    <INPUT TYPE="hidden" NAME="request" VALUE=
    and the word after VALUE is your request

To find out the WHATIF output name for the server that you need: 
  - Go to http://swift.cmbi.kun.nl/WIWWWI/ and click on your server 
  - Submit a single job
  - Click on Results
  - Look for an hyperlinked (blue-underlined) filename, e.g. fixed.pdb or hadded.pdb
        
"""
import spyder

site = 'http://swift.cmbi.ru.nl/'
whatifcgi = 'wiw-cgi/FiledownCGI.py'

"""
  THESE ARE THE LOCATIONS OF THE WHATIF SITE AND SERVER, THIS MAY CHANGE IN THE FUTURE!!
"""


import os, sys, time

if spyder.python2:
  import imp
  http = imp.new_module("http")
  import httplib, urllib, urlparse
  http.client = httplib
  urllib.parse = urlparse
  import mimetypes
else:
  import http.client, urllib, urllib.parse
  import mimetypes
   
def post_multipart(filevalue, request):
    spyder.load("http")
    url = site + whatifcgi
    fields = (('request', request), ('&PDB1', ''))
    filekey = '&FIL1'

    urlparts = urllib.parse.urlsplit(url)
    host = urlparts[1]
    selector = urlparts[2]
    
    fil = (filekey, "dummy", filevalue)
    
    content_type, body = spyder.http.httppost.encode_multipart_formdata(fields, [fil])
    counter = 0
    while 1:
      counter += 1
      ok = True
      try:
        h = http.client.HTTPConnection(host)
        headers = {
          'User-Agent': 'anonymous',
          'Content-Type': content_type    
        }
        h.request('POST', selector, body, headers)
        res = h.getresponse() 
      except:
        ok = False
        if counter == 5: raise 
      if ok: break
    return res.read()

    
def whatif(pdbdata, request, outputname): 
  url = site + whatifcgi
  
  
  urlparts = urllib.parse.urlsplit(url)
  host = urlparts[1]
  selector = urlparts[2]
    
  try:
    loc = post_multipart(pdbdata, request)
    fields = []
    for l in loc.splitlines():
      l = l.strip()
      if l.startswith("<INPUT TYPE"):
        name = l[l.index("NAME=")+len("NAME="):].strip()[1:]
        name = name[:name.index('"')]
        val = l[l.index("VALUE=")+len("VALUE="):].strip()[1:]
        val = val[:val.index('"')]
        fields.append((name,val))
    for f in fields:
      if f[0] == 'ID':
        location = "/servers/tmp/" + f[1]
        break
    
  except:
    raise spyder.atom.AtomValidationError("Could not submit job to the WHATIF server\n")
  
  try:  
    return urllib.urlopen(site+ location + "/" + outputname).read()
  except:
    raise spyder.atom.AtomValidationError("Cannot download WHATIF file")
  
def add_hydrogens(pdbdata):
  return whatif(pdbdata,"htopo", "hadded.pdb")

def analyze_protonation_state(pdbdata,pdbname="the PDB file"):      
  class protonationstate:
    pass
  lines = pdbdata.splitlines()
    
  histidines = {}
  
  for l in lines:
    code = l[17:20]
    if code == "HIS":
      nr = int(l[22:26])
      if nr not in histidines:
        p = protonationstate()
        p.d1 = False
        p.d2 = False      
        p.e1 = False
        p.e2 = False
        histidines[nr] = p
  
  hisprotonatoms = (" HD1"," HD2"," HE1"," HE2")      
        
  for l in lines:
    code = l[17:20]
    atom = l[12:16]
    if code == "HIS" and atom in hisprotonatoms:
      nr = int(l[22:26])
      if nr not in histidines: continue
      currhis = histidines[nr]
      if atom == hisprotonatoms[0]: currhis.d1 = True
      if atom == hisprotonatoms[1]: currhis.d2 = True
      if atom == hisprotonatoms[2]: currhis.e1 = True
      if atom == hisprotonatoms[3]: currhis.e2 = True
  
  ret = []
  for nr in histidines:
    his = histidines[nr]
    dcount = his.d1 + his.d2
    ecount = his.e1 + his.e2
    totcount = dcount + ecount
    if totcount == 3:
      if dcount == 2:       
        ret.append(dict(resid=nr,state="HISD"))
      else: ret.append(dict(resid=nr,state="HISE"))
    elif totcount == 4:
      ret.append(dict(resid=nr,state="HIS+"))
    else: 
      raise spyder.atom.AtomValidationError("WHATIF could not guess the protonation state of histidine %d in %s" % (nr, pdbname))
  return ret
      
    
