import sys, os
topdir = ".." + os.sep + os.path.split(__file__)[0]
sys.path.append(topdir)
import attracthtmlform 
import spyder
from attractmodel import AttractModel
import form

cgi = "http://localhost/cgi/something.py"
if len(sys.argv) > 1: cgi = sys.argv[1]

f = AttractModel._form()
f = form.webform(f, partnerlength=1) #generate HTML for only one partner, the server .js will clone
html = attracthtmlform.htmlform(form=f, cgi=cgi, header=form.header, footer=form.footer, header_indentation = 12)
print(html)