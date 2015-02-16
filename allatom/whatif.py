"""
Python script to contact the WHATIF web servers at the CMBI
In principle, this will work with any WHATIF server that accepts a single file and outputs a single file

Sjoerd de Vries, 2006, 2007, 2014
          
To find out the WHATIF request for the server that you need: 
  - Go to http://swift.cmbi.ru.nl/WIWWWI/ and click on your server 
  - Click on View Frame Source
  - Look for this line:
    <INPUT TYPE="hidden" NAME="request" VALUE=
    and the word after VALUE is your request

To find out the WHATIF output name for the server that you need: 
  - Go to http://swift.cmbi.ru.nl/WIWWWI/ and click on your server 
  - Submit a single job
  - Click on Results
  - Look for an hyperlinked (blue-underlined) filename, e.g. fixed.pdb or hadded.pdb
        
"""

from __future__ import print_function

site = 'http://swift.cmbi.ru.nl/'
whatifcgi = 'wiw-cgi/FiledownCGI.py'
whatifcgi2 = 'wiw-cgi/ChecksCGI.py'
whatifcgi3 = 'wiw-cgi/GenericCGI.py'

"""
  THESE ARE THE LOCATIONS OF THE WHATIF SITE AND SERVER, THIS MAY CHANGE IN THE FUTURE!!
"""


import os, sys, time

python2 = (sys.version_info[0] == 2)
python3 = (sys.version_info[0] == 3)

if python2:
  import imp
  http = imp.new_module("http")
  import httplib, urllib, urlparse
  http.client = httplib
  urllib.parse = urlparse
  import mimetypes
else:
  import http.client, urllib, urllib.parse
  import mimetypes

def post_multipart(url, fields, files):
    urlparts = urllib.parse.urlsplit(url)
    host = urlparts[1]
    selector = urlparts[2]
        
    content_type, body = encode_multipart_formdata(fields, files)
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
  
def encode_multipart_formdata(fields, files):
    """
    Based on some code snippet found on the web in 2006 or so
    """
    BOUNDARY = '----------1234567890abcdefghij_$'
    CRLF = '\r\n'
    L = []
    for file in files:
      (filekey, filename, filevalue) = file
      L.append('--' + BOUNDARY)
      L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (filekey, filename))
      L.append('Content-Type: text/plain')
      L.append('')
      for line in filevalue.split('\n'):
        L.append(line)    
    for (key, value) in fields:
        if value == None: value = ""
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % str(key))
        L.append('')
        L.append(str(value))
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body

###############################################################################
    
def whatif0(cgi, pdbdata, request, outputname): 

  url = site + cgi
  fields0 = (('request', request), ('&PDB1', ''))
  filekey = '&FIL1'  
  fil = (filekey, "dummy", pdbdata)
    
  try:
    loc = post_multipart(url, fields0, [fil])
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
    raise Exception("Could not submit job to the WHATIF server\n")
  
  try:  
    return urllib.urlopen(site+ location + "/" + outputname).read()
  except:
    raise Exception("Cannot download WHATIF file")

def whatif(pdbdata, request, outputname): 
  return whatif0(whatifcgi, pdbdata, request, outputname)

def whatif2(pdbdata, request, outputname): 
  return whatif0(whatifcgi2, pdbdata, request, outputname)
  
def whatif3(pdbdata, request):
  """
  For WHATIF services that do not return a file, but direct results
  """
  url = site + whatifcgi
  fields0 = (('request', request), ('&PDB1', ''))
  filekey = '&FIL1'  
  fil = (filekey, "dummy", pdbdata)
  
  try:
    job = post_multipart(url, fields0, [fil])
    fields = []
    for l in job.splitlines():
      l = l.strip()
      if l.startswith("<INPUT TYPE"):
        name = l[l.index("NAME=")+len("NAME="):].strip()[1:]
        name = name[:name.index('"')]
        val = l[l.index("VALUE=")+len("VALUE="):].strip()[1:]
        val = val[:val.index('"')]
        fields.append((name,val))
    location = None
    for f in fields:
      if f[0] == 'ID':
        location = f[1]
        break
    assert location is not None    
  except:
    raise Exception("Could not submit job to the WHATIF server\n")
  
  url2 = site + whatifcgi3
  fields = [("ID", f[1]), ("refresh", request)]
  result = post_multipart(url2, fields, [])
  return(result)
  
###############################################################################  

services = {
  "htopo": ("Add hydrogens", "whatif", "hadded.pdb"),
  "corall": ("Complete sidechains", "whatif", "fixed.pdb"),
  "plnchk": ("Check planarity", "whatif3"),
  "ramchk": ("Ramachandran plot evaluation", "whatif3"),
  "angchk": ("Anomalous bond angles", "whatif3"),
  "bndchk": ("Anomalous bond lengths", "whatif3"),
  "omechk": ("Omega angle distribution", "whatif3"),
  "modcheck": ("Protein model check (human-readable)", "whatif2", "pdbout.txt"),
  "modcheck2": ("Protein model check (computer-parsable)", "whatif2", "check.db"),
}
  
def run(service, pdbdata): 
  s = services[service]
  command = s[1]
  assert command in ("whatif", "whatif2", "whatif3"), command
  if service == "modcheck2": service = "modcheck"
  if command in ("whatif", "whatif2"): 
    assert len(s) == 3, s
    result = whatif(pdbdata, service, s[2])
  else:
    assert len(s) == 2, s
    result = whatif2(pdbdata, service)
  return result

if __name__ == "__main__":
  import sys
  try:
    service = sys.argv[1]
    pdbfile = sys.argv[2]
    pdbdata = open(pdbfile).read()
    s = services[service]
  except:
    print("*" * 50)
    print("ERROR: ")
    print("*" * 50)
    import traceback
    traceback.print_exc()
    print("*" * 50)
    print("Usage: whatif.py <service> <pdb file>")
    print("*" * 50)
    print("Possible services:")
    for k,v in sorted(services.items()):
      print(k + " : " + v[0])
    print("*" * 50)
    sys.exit()
  result = run(service, pdbdata)  
  print(result)  
  