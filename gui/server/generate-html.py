import sys, os
import attracthtmlform 
import spyder
from attractmodel import AttractModel
import form

cgi = sys.argv[1]

f = AttractModel._form()
f = form.webform(f)
html = attracthtmlform.htmlform(
 form=f, cgi=cgi, 
 header=form.header, footer=form.footer, header_indentation = 12, 
 newtab=True
)
print(html)