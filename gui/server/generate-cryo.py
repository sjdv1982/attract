import sys, os
sys.path.insert(0, "..")
import spyder
import attractmodel
from cryo import CryoPartnerRun
import form_cryo

cgi = sys.argv[1]

f = CryoPartnerRun._form()
f = form_cryo.webform(f)
html = form_cryo.html(f, cgi, None, newtab=True)
print(html)