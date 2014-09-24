import sys, os
sys.path.insert(0, "..")
import spyder
spyder.silent = True
import attractmodel
from attractmodel import NARefineInterface
import form_narefine

cgi = sys.argv[1]

f = NARefineInterface._form()
f = form_narefine.webform(f)
html = form_narefine.html(f, cgi, None, newtab=True)
print(html)