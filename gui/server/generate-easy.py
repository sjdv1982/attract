import sys, os
sys.path.insert(0, "..")
import spyder
from attractmodel import AttractEasyModel
import formeasy, form_model

cgi = sys.argv[1]

f = AttractEasyModel._form()
f = formeasy.webform(f)
html = form_model.html(f, cgi, newtab=True)
print(html)