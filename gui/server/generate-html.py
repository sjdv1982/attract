import sys, os
import spyder
from attractmodel import AttractModel
import form_model

cgi = sys.argv[1]

f = AttractModel._form()
f = form_model.webform(f)
html = form_model.html(f, cgi, newtab=True)
print(html)