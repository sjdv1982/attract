from ..Handler import Constructor, CMethod, XMethod, Method, UniMethod, PassDown, NullHandler
from . import tag

class Printer(object):
  def __init__(self):
    self.clear()
  def clear(self):
    self.value = ""
    self.indentation = 0
    self._has_indented = False
  def write(self, s):
    if not self._has_indented:
      self.value += " " * self.indentation
      self._has_indented = True
    self.value += s
  def newline(self):
    self.value += "\n"
    self._has_indented = False
  def indent(self, increment):
    self.indentation += increment
  def dedent(self, decrement):
    self.indentation -= decrement

def init(root):
  #normal markup tags may only contain markup children (markup tag and string tag)
  root.addHandler( Constructor("markup-tag", tag.init_tag) )
  root.addHandler( NullHandler("markup-tag", "make-markup") )
  root.addHandler( NullHandler("markup-tag", "complete-markup") )
  root.addHandler( XMethod("markup-tag", "print-markup", tag.print_tag) )
  
  #mixed tags may not contain non-markup children
  #complete-markup will convert it into a normal markup tag
  root.addHandler( Constructor("markup-mixedtag", tag.init_tag) )
  root.addHandler( PassDown("markup-mixedtag", "make-markup") )
  root.addHandler( XMethod("markup-mixedtag", "complete-markup", tag.mixedtag_complete_markup) )
  
  #string "tag" (text content inside a tag)
  root.addHandler( Constructor("markup-stringtag", tag.init_stringtag) )
  root.addHandler( NullHandler("markup-stringtag", "make-markup") )
  root.addHandler( NullHandler("markup-stringtag", "complete-markup") )
  root.addHandler( Method("markup-stringtag", "print-markup", tag.print_stringtag) )

  #null tag (will be replaced when markup is completed)
  root.addHandler( Constructor("markup-nulltag", tag.init_nulltag) )
  root.addHandler( PassDown("markup-nulltag", "make-markup") )
  root.addHandler( XMethod("markup-nulltag", "complete-markup", tag.nulltag_complete_markup) )
  
  root.space.printer = Printer()
  