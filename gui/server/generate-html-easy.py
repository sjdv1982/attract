import sys, os
sys.path.insert(0, "..")
import spyder
spyder.silent = True
from attractmodel import AttractEasyModel
import formeasy

cgi = sys.argv[1]

f = AttractEasyModel._form()
f = formeasy.webform(f)
html = formeasy.html(f, cgi, None, newtab=True)
print(html)