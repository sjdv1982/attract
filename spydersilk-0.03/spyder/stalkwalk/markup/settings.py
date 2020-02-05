from ..Space import Space

class BaseSettings(Space):
  _slots = () #emulates __slots__
  def __init__(self):
    Space.__init__(self, top=False)
  def __setattr__(self, attr, value):
    if attr not in self._slots: raise AttributeError(attr)
    Space.__setattr__(self, attr, value)
  def validate(self):
    pass
  def update(self, d):
    for key in d: setattr(self, key, d[key])

class AttributeSettings(BaseSettings):
  _slots = (
    "quotes",
    "minimized",
    "space_before",
    "space_after",
    "space_before_assign",
    "space_after_assign",    
  )
  def validate(self):
    assert self.quotes in (True, False, None), self.quotes
    assert self.minimized in (True, False, None), self.minimized
    assert not self.quotes or not self.minimized    
    assert self.space_before is None or isinstance(self.space_before, int), self.space_before
    assert self.space_after is None or isinstance(self.space_after, int), self.space_after
    assert self.space_before_assign is None or isinstance(self.space_before_assign, int), self.space_before_assign
    assert self.space_after_assign is None or isinstance(self.space_after_assign, int), self.space_after_assign
    if self.minimized: assert not self.space_before_assign and not self.space_after_assign
    
class TagSettings(BaseSettings):    
  _slots = (
    "inline",    
    "indent",
    "closetag",
    "lines_before",
    "lines_after",
    "lines_inside_before",
    "lines_inside_after",    
    "comment",
    "anti", #reverse opentag and closetag
  )
  def validate(self):
    assert self.inline in (True, False, None), self.inline
    assert self.indent in (True, False, None), self.indent
    assert self.indent is None or self.inline != True, (self.indent, self.inline)
    assert self.closetag in (True, False, "slash", None), self.closetag
    if self.inline: assert self.closetag != False
    assert self.lines_before is None or isinstance(self.lines_before, int), self.lines_before
    assert self.lines_after is None or isinstance(self.lines_after, int), self.lines_after
    assert self.comment is None or isinstance(self.comment, str), self.comment
    if self.comment is not None:
      assert self.closetag in (False, "slash"), self.closetag
      assert self.comment.lstrip().startswith("<!--"), self.comment
      assert self.comment.rstrip().startswith("-->"), self.comment
    assert self.anti in (True, False, None), self.anti
    if self.anti:
      assert self.closetag in (True, None)
    assert self.lines_inside_before is None or isinstance(self.lines_inside_before, int), self.lines_inside_before
    assert self.lines_inside_after is None or isinstance(self.lines_inside_after, int), self.lines_inside_after
    if self.lines_inside_before or self.lines_inside_after: 
      assert self.inline == False
      assert self.closetag in (None, True)

class OCTagSettings(BaseSettings):
  _slots = (
    "space_before", #padding between < and tag name
    "space_after",#padding between tag name or last attribute and > or />
    "comment",
  )
  def validate(self):
    assert self.comment is None or isinstance(self.comment, str), self.comment
    assert self.space_before is None or isinstance(self.space_before, int), self.space_before
    assert self.space_after is None or isinstance(self.space_after, int), self.space_after
    if self.comment is not None:
      assert self.comment.lstrip().startswith("<!--"), self.comment
      assert self.comment.rstrip().startswith("-->"), self.comment
  
    