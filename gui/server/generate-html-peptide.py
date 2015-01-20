import sys, os
sys.path.insert(0, "..")
import spyder
spyder.silent = True
from attractmodel import AttractPeptideModel
import form_peptide

cgi = sys.argv[1]

f = AttractPeptideModel._form()
f = form_peptide.webform(f)
html = form_peptide.html(f, cgi, None, newtab=True)
print(html)