from .settings import AttributeSettings, TagSettings, OCTagSettings
from .. import Handler

def validate_params(params):
  """
  Params validation for both markuptag and mixedtag
  """
  assert params.tag is not None
  assert isinstance(params.tag, tuple) and len(params.tag) == 2, params.tag
  assert isinstance(params.tag[0], str), params.tag
  assert isinstance(params.tag[1], TagSettings), params.tag
  
  assert params.opentag is not None
  assert isinstance(params.opentag, OCTagSettings), params.opentag
  
  if params.closetag is not None:
    assert params.tag[1].closetag in (None, True), params.tag[1].closetag
    assert isinstance(params.closetag, OCTagSettings), params.closetag
  else:
    assert params.tag[1].closetag in ("slash", False), params.tag[1].closetag
  
  if params.attributes is not None:
    for attribute in params.attributes:
      assert isinstance(attribute, tuple) and len(attribute) == 3, attribute
      assert isinstance(attribute[0], str), attribute
      assert isinstance(attribute[2], AttributeSettings), attribute
  
def init_tag(node):
  """
  Constructor for both markup tags and mixed tags
  """
  space = node.space
  params = node.parent.space.params
  
  space.tag = params.tag
  space.opentag = params.opentag
  space.closetag = params.closetag
  space.attributes = params.attributes
   
  validate_params(params)

def print_octag(printer, tagname, tag, 
  opentag, slashclose = False, attributes = []
 ):  
  #prints open/close tag
  s = "<"
  if not opentag: s += "/"
  if tag.space_before: s += tag.space_before * " "
  s += tagname  
  
  if attributes:
    for aname, avalue, asettings in attributes:
      if avalue is None and not asettings.minimized: continue
      s += " "
      if asettings.space_before: s += asettings.space_before * " "
      s += aname
      if not asettings.minimized:
        if asettings.space_before_assign: 
          s += asettings.space_before_assign * " "
        s += "="
        if asettings.space_after_assign: 
          s += asettings.space_after_assign * " "            
        if asettings.quotes: s += '"'
        s += str(avalue)
        if asettings.quotes: s += '"'
      if asettings.space_after: s += asettings.space_after * " "  
      

  if tag.space_after: s += tag.space_after * " "
  if slashclose: s += "/"
  s += ">"
  if tag.comment: s += tag.comment
  printer.write(s)  
    
def print_tag(message, node):
  space = node.space
  params = node.parent.space.params
  printer = params.printer
  if printer is None: printer = node.root.space.printer
  assert printer is not None
  
  tagname = space.tag[0]
  indent = space.tag[1].indent
  if indent is None: indent = node.root.space.indent
  if indent is None: indent = 2 
  inline = space.tag[1].inline
  if inline is None: inline = False
  anti = space.tag[1].anti
  if anti is None: anti = False
    
  if space.tag[1].lines_before:
    for n in range(space.tag[1].lines_before):
      printer.newline()
  
  slashclose = (space.tag[1].closetag == "slash")
  
  if anti:
    print_octag(printer, tagname, space.closetag, False)
  else:
    print_octag(printer, tagname, space.opentag, True, slashclose, space.attributes)
  
  if space.closetag is not None and not inline:
    printer.newline()
    if space.tag[1].lines_inside_before:
      for n in range(space.tag[1].lines_inside_before):
        printer.newline()
    printer.indent(indent)      
  
  space.params.printer = printer
  for childname in list(node.getChildren()):
    if node.hasChild(childname): #could be removed now...
      child = node[childname]
      child.sendMessage(message)
      if not space.tag[1].inline: printer.newline()
  space.params.clear()
  
  if space.closetag is not None and not inline:
    if space.tag[1].lines_inside_after:
      for n in range(space.tag[1].lines_inside_after):
        printer.newline()    
    printer.dedent(indent)
  
  if anti:
    print_octag(printer, tagname, space.opentag, True, False, space.attributes)
  else:
    if space.closetag is not None:
      print_octag(printer, tagname, space.closetag, False)
    
  if space.tag[1].comment:
    printer.write(space.tag[1].comment)    
  if space.tag[1].lines_after:
    for n in range(space.tag[1].lines_after):
      printer.newline()
      
def mixedtag_complete_markup(message, node):
  Handler.passdown(message, node)  
  node.morph("markup-tag")
  
def init_stringtag(node):  
  space = node.space
  params = node.parent.space.params
  
  space.value = params.value
  assert params.value is not None, (node.parent.name, node.parent.space.formpath)
  assert isinstance(params.value, str), params.value
    
def print_stringtag(node):  
  space = node.space
  params = node.parent.space.params
  
  params.printer.write(space.value)
  
def init_nulltag(node): pass
def nulltag_complete_markup(message, node):
  Handler.passdown(message, node)  
  ch = node.getChildren()
  if len(ch) == 0:
    node.parent.removeChild(node.name)
  elif len(ch) == 1:
    child = node.detachChild(node.getChildren()[0])
    node.parent.replaceChild(node.name, child)
  else:
    p = node.parent
    mypos = node.parent.getChildren().index(node.name)
    children = list(node.getChildren())
    node.parent.detachChild(node.name)
    for childnr,childname in enumerate(children):
      child = node.detachChild(childname)
      p.attachChild(child, str(childname)+"-"+str(mypos), mypos+childnr)
    node.sendMessage("destroy")
    node.sendMessage("post-destroy")