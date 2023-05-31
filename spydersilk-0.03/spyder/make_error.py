from __future__ import print_function
import sys, os
import spyder, Spyder
from spyder.spydercompile import getlines, getblock, quotematch
from spyder.modules.core import parse_error
        
spydertypename = sys.argv[1]
if len(sys.argv) > 2:
  for a in sys.argv[2:]:
    if a.find(os.sep) > -1:
      sys.path.insert(0, a)
    else:
      __import__(a)    
    
spyder.load(spydertypename)
spydertype = getattr(Spyder, spydertypename)

s = ""
doclines = spydertype.__doc__.splitlines()
mode = False
lnr = -1
while 1:
  lnr += 1
  if lnr == len(doclines): break
  l = doclines[lnr]
  sepline = (l.find("-"*20) > -1)    
  if mode:
    if sepline: break
    s += l + "\n"
  else:
    if sepline: 
      mode = True
      lnr += 1
    continue

error = spyder.core.error[spydertypename] 
done = set()
lines = getlines(s)
for l in lines:
  name,title,block,blockcomment = getblock(l)
  if name != "Type": continue
  assert title == spydertypename or title.startswith(spydertypename+"("), title
  tblocks = getlines(block)
  for tblock in tblocks:
    sname,stitle,sblock,sblockcomment = getblock(tblock)
    if sname == "error":
      #print("ERROR BLOCK", sblock)
      block_errors = parse_error(sblock)
      for statement,newstatement in block_errors:
        if statement not in error:
          print("WARNING: unknown error message: '%s'" % statement, file=sys.stderr)
          continue          
        assert statement in error, statement
        done.add(statement)
  
errorblock = []
for statement in error:
  if statement not in done: errorblock.append(statement)    

slines = []
slines.append("  error {")
for statement in errorblock:
  p = "    "
  p += repr(statement) + '\n'
  p += "    =>\n"
  p += '    \'\'\n'
  slines.append(p)
slines.append("  }")
print("\n".join(slines))
