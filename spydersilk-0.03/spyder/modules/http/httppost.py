# Copyright 2008, 2009 Sjoerd de Vries
# This file is part of the Spyder module: "http" 
# For licensing information, see LICENSE.txt 

"""
Library for remote web access (early development stage)
"""

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
