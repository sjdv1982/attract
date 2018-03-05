import sys, os
sys.path.insert(0, "..")
import spyder
spyder.silent = True
import attractmodel
from attractmodel import CryoEasyInterface
import form_cryo_easy

cgi = sys.argv[1]

f = CryoEasyInterface._form()
f = form_cryo_easy.webform(f)
html = form_cryo_easy.html(f, cgi, None, newtab=True)
print(html)
