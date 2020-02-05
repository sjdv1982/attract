from .. import Handler, Space
from ..hamlet.tag import SimpleTag, VoidTag, ComplexTag, ComplexVoidTag, TextWrapTag, StringTag, NullTag


def _add_headers(node, tag, begin=False):
  space = node.space
  if space.othertokens is None: return
  headers = [h[1][0] for h in space.othertokens if h[0] == "header"]  
  for hnr,h in enumerate(headers):    
    header = ComplexTag(
     "tr", 
     children = [ TextWrapTag("td", h, [("class", "h2"), ("colspan", 2)], inline = True) ],
     inline=True,
    )
    if begin:
      tag.attachChild(header, "header-"+str(hnr+1), hnr)
    else:
      tag.attachChild(header, "header-"+str(hnr+1))
      
def _container_tag(node, tag):
  assert node.space.htmlnode is None, (node.space.formpath, node.space.htmlnode)
  assert tag is not None, node.space.formpath
  Handler.passdown("html-make-concrete", node)
  node.space.htmlnode = tag  
  for childname in node.getChildren():
    child = node[childname]
    if child.space.complexity is not None:
      tip = child.space.htmltip
      if tip is None: tip = child.space.htmlnode      
      assert tip is not None, child.space.formpath
      _add_headers(child, tip, begin=True)
      continue
    nulltag = NullTag()
    _add_headers(child, nulltag)
    h = str(child.space.formname)
    head = TextWrapTag("td", h, inline = True)    
    trtag = SimpleTag("tr", lines_after=1)
    tdtag = SimpleTag("td", [("class", "formdata")])
    trtag.attachChild(head, "head")
    trtag.attachChild(tdtag, "child")
    trtag.space.htmltip = tdtag
    nulltag.attachChild(trtag, "child")
    node.detachAndReplaceChild(childname, nulltag)        
    tdtag.attachChild(child, "child")
  
  node.morph("html-concrete-node")
  return Handler.CLAIM

def morph_formwrapper(node):
  assert node.space.cgi is not None
  attributes = [
    ("encType", "multipart/form-data"),
    ("action", node.space.cgi),
    ("method", "post"),
  ]
  if node.space.newtab:
    attributes.append(("target", "_blank"))  
  formtag = SimpleTag("form", attributes, inline=True)
  formtag.attachChild(node.space.htmlnode, "table")
  node.space.htmlnode = None
  _container_tag(node, formtag)
  
def morph_h1wrapper(node):
  formname = node.space.formname
  if node.space.formnode is not None:
    formname = node.space.group.name
  nulltag = NullTag()
  n = NullTag()
  trtag = ComplexTag (
   "tr",
   children = [ 
    TextWrapTag("td", formname, [("class","h1"),("colspan",2)], inline = True) 
   ],
   inline = True
  )
  nulltag.attachChild(trtag, "child1")
  nulltag.attachChild(n, "child2")
  nulltag.space.htmltip = n
  return _container_tag(node, nulltag)

def morph_foldwrapper(node):
  formpath = node.space.formpath
  formname = node.space.formname
  if node.space.formnode is not None:
    formpath = node.space.formnode.space.formpath + ("group" + str(node.space.groupindex+1),)
    formname = node.space.group.name
  assert formpath is not None, (node.name, node.space.formpath, node.nodetype)
  fp = "-".join([str(s) for s in formpath[1:]])  
  tabletag = SimpleTag("table", anti=True, inline=True)
  divtag = SimpleTag(
   "div", 
   [("class","topItem"), ("id","level%d" % node.space.wrapping )], 
   inline = True 
  )
  tabletag.attachChild(divtag, "child")
  lines_before = 0
  if node.space.wrapping == 1 or node.space.formnode is None: lines_before = 1
  if node.space.wrapping == 1 and node.space.formnode is not None: lines_before = 0
  spantag = ComplexTag("span", 
   [
     ("id","foldmenu-%s-title" % fp), 
     ("style","float:right"), 
   ],     
   [VoidTag("img",(("src","arrow-down.png",{"space_after":1}),), closetag = "slash")],
   inline = True,
   lines_before = lines_before,
  )
  divtag.attachChild(spantag, "child")  
  stag = StringTag(formname)
  divtag.attachChild(stag, "child2")
  div2tag = ComplexTag (
    "div",
    [
      ("id","foldmenu-%s" % fp),
      ("class","switchgroup1"),
    ],
    [ SimpleTag("table") ],
    inline = True
  )
  tabletag.attachChild(div2tag, "child2")
  tabletag.space.htmltip = div2tag 
  return _container_tag(node, tabletag)

def morph_divwrapper(node):
  divtag = SimpleTag("div")
  return _container_tag(node, divtag)
  
def _leaf_tag(node, tag):
  node.space.htmlnode = tag
  Handler.passdown("html-make-concrete", node)  
  node.morph("html-concrete-node")
  return Handler.CLAIM

def _extract_fields(space, form, variables):
  ret = Space.Space(top=False)
  for v in variables:
    if getattr(space, v) is not None:
      setattr(ret, v, getattr(space, v))
    elif form is not None and getattr(form, v) is not None:
      setattr(ret, v, getattr(form, v))
  return ret    

def morph_inputtag_text_spin(node):
  space = node.space
  if space.formtype != "spin": return Handler.NO_MATCH
  return _morph_inputtag_text(node, "text")

def morph_inputtag_text(node):
  space = node.space
  if space.formtype not in ("text", "password"): return Handler.NO_MATCH
  return _morph_inputtag_text(node, space.formtype)
  
def morph_textarea(node):
  space = node.space
  if space.formtype != "textarea": return Handler.NO_MATCH
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  fields = ("defaultvalue", "rows")
  p = _extract_fields(space, space.form, fields)
  if p.defaultvalue is None: p.defaultvalue = ""
  p.defaultvalue = str(p.defaultvalue) 
  attributes = [
    ("name", fp),
    ("rows", p.rows),
    ("id", fp),
  ]  
  textareatag = TextWrapTag("textarea", p.defaultvalue, attributes, inline=True)
  return _leaf_tag(node, textareatag)

def _morph_inputtag_text(node, typ):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  fields = ("defaultvalue", "length", "placeholder")
  p = _extract_fields(space, space.form, fields)
  if p.defaultvalue is None: p.defaultvalue = ""
  p.defaultvalue = str(p.defaultvalue) 
  attributes = [
    ("type", typ), 
    ("name", fp),
    ("id", fp),
    ("value", p.defaultvalue, {"space_before_assign": 1, "space_after_assign": 1} ),
    ("maxlength", p.length),
    ("placeholder", p.placeholder),
  ]  
  inputtag = VoidTag("input", attributes)
  return _leaf_tag(node, inputtag)

def morph_selecttag(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "option": return Handler.NO_MATCH  
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("name", fp),
    ("id", fp),
  ]  
  options, optiontitles = space.options, space.optiontitles
  assert len(options) and len(options) == len(optiontitles)
  children = []
  if p.defaultvalue is None:
    child = SimpleTag("option", (("value",""),"selected"), inline=True )
    children.append(child)
  for opt, optitle in zip(options, optiontitles):  
    attrs = (("value",opt),)
    if opt == p.defaultvalue: attrs = (("value",opt),"selected") 
    child = TextWrapTag("option", optitle, attrs, inline=True)
    children.append(child)
  inputtag = ComplexTag("select", attributes, children)
  
  return _leaf_tag(node, inputtag)

def morph_inputtag_file(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "file": return Handler.NO_MATCH
  fields = ("placeholder",)
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("type", "file"), 
    ("name", fp),
    ("id", fp),
    ("value", "", {"space_before_assign": 1, "space_after_assign": 1}),
    ("placeholder", p.placeholder),
  ]  
  inputtag = VoidTag("input", attributes)
  return _leaf_tag(node, inputtag)

def morph_inputtag_checkbox(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "checkbox": return Handler.NO_MATCH
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("type", "checkbox"), 
    ("name", fp),
    ("id", fp),
    ("value", "", {"space_before_assign": 1, "space_after_assign": 1}),
  ]
  if p.defaultvalue: attributes.append("checked")
  inputtag = VoidTag("input", attributes)
  return _leaf_tag(node, inputtag)

def morph_inputtag_hidden(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "hidden": return Handler.NO_MATCH
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  if p.defaultvalue:
    raise ValueError("Data payload in a hidden HTML field not yet supported: %s" % str(node.space.formpath))
  node.parent.detachChild(node.name)
  return Handler.CLAIM
    
def morph_inputtag_unknown(node):
  raise ValueError("Unknown form type", (node.space.formtype, node.space.formpath))