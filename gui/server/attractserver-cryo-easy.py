#!/usr/bin/env python2
from __future__ import print_function

def report(s=""):
  print(str(s).replace("[ADVANCED]","").replace("[/ADVANCED]",""))

report("Content-type: text/html\n\n<pre>")

try:
  import os, sys, time, traceback
  sys.path.append(os.path.split(os.path.abspath(__file__))[0] + "/..")
  import spyder, Spyder
  import spyder.formtools
  import attractmodel
  import form_cryo_easy
  import attractsave
  from serverlib import serve_attract, AttractServerError
  os.system("chmod a+r+w /tmp/*spy.*")
except Exception, e:
  s = "<b>There was an error in the initialization of the server</b>\n"
  s += """Please save this page from your browser as HTML and email it to sjoerd@tum.de, and it will be fixed as soon as possible.
The ATTRACT web interface is in active development, thank you for your patience.
"""
  s += "\n<B>Full error information</B>:\n\n"
  s += traceback.format_exc()
  report(s)
  report("</pre>")
  sys.exit()

try:
  r = serve_attract(Spyder.CryoEasyInterface, form_cryo_easy, attractsave.deploy_cryo_easy, easy=False)
  report(r)
except AttractServerError as e:
  report("<b>There was an inconsistency in your data</b>")
  report("\n<B>Error message</B>")
  report(e.status)
  report()
  if e.delta is not None: report(e.delta)
  report("</pre>")
  sys.exit()
except Exception as e:
  report("There was an unknown error in the server")
  report("""Please save this page from your browser as HTML and email it to sjoerd@tum.de, and it will be fixed as soon as possible.
The ATTRACT web interface is in active development, thank you for your patience.
""")
  report(traceback.format_exc() + "\n")
  report("</pre>")
  sys.exit()

report("</pre>")
