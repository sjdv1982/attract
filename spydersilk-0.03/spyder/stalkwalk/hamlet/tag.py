from ..markup.settings import TagSettings, OCTagSettings, AttributeSettings
from . import node as hnode
from ..Node import Node

def parse_attributes(attributes):
  if attributes is None: return None
  ret = []
  for at in attributes:
    iterable = True
    settings = AttributeSettings()
    try:
      len(at)
    except TypeError:
      iterable = False  
    
    if isinstance(at, str):      
      settings.minimized = True
      at2 = (at, None, settings)
    elif iterable and len(at) == 2:
      if isinstance(at[1], AttributeSettings):
        assert at[1].minimized, at
        at2 = (at[0], None, at[1])
      elif isinstance(at[1], dict):
        assert at[1].minimized, at
        settings.update(at[1])
        at2 = (at[0], None, settings)        
      else:
        if isinstance(at[1], str):
          settings.quotes = True
        at2 = (at[0], at[1], settings)
    elif iterable and len(at) == 3:      
      if isinstance(at[2], AttributeSettings):
        at2 = at
      else: #dict or dict-like
       if isinstance(at[1], str):
         settings.quotes = True      
       settings.update(at[2])
       at2 = (at[0], at[1], settings)
    else:
      raise TypeError("Invalid attribute format", at)
    ret.append(at2)
  ret2 = []
  for a1,a2,a3 in ret:
    if a1.endswith("_"): a1 = a1[:-1]
    ret2.append((a1,a2,a3))
  return ret2

def _tag(hnode):  
  node = Node("markup-mixedtag", None, hnode, None)
  node.sendMessage("init")
  hnode.space.params.clear()
  node.parent = None
  return node
  
def SimpleTag(tag, attributes=None, **kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.closetag = OCTagSettings()
  params.attributes = parse_attributes(attributes)
  return _tag(hnode)

def VoidTag(tag, attributes=None, **kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.closetag = False  
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.attributes = parse_attributes(attributes)
  return _tag(hnode)

def StringTag(text):  
  params = hnode.space.params
  params.value = text
  node = Node("markup-stringtag", None, hnode, None)
  node.sendMessage("init")
  hnode.space.params.clear()
  node.parent = None
  return node
  
def TextWrapTag(tag, text, attributes=None, **kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.closetag = OCTagSettings()
  params.attributes = parse_attributes(attributes)
  tag = _tag(hnode)
  tag.space.params.value = text
  tag.addChild("markup-stringtag", "text")
  tag.space.params.clear()  
  return tag

def ComplexTag(tag, attributes=None, children=[],**kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.closetag = OCTagSettings()
  params.attributes = parse_attributes(attributes)
  tag = _tag(hnode)
  for childnr, child in enumerate(children):
    assert isinstance(child, Node)
    tag.attachChild(child, str(childnr), childnr)  
  return tag

def ComplexVoidTag(tag, attributes=None, children=[],**kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.attributes = parse_attributes(attributes)
  tag = _tag(hnode)
  for childnr, child in enumerate(children):
    assert isinstance(child, Node)
    tag.attachChild(child, str(childnr), childnr)  
  return tag
  
def ContainerTag(tag, children=[],**kwargs):  
  params = hnode.space.params
  settings = TagSettings()
  settings.update(kwargs)
  params.tag = (tag, settings)
  params.opentag = OCTagSettings()
  params.closetag = OCTagSettings()
  tag = _tag(hnode)
  for childnr, child in enumerate(children):
    assert isinstance(child, Node)
    tag.attachChild(child, str(childnr), childnr)  
  return tag
  
def NullTag(children = []):    
  node = Node("markup-nulltag", None, hnode, None)
  node.sendMessage("init")
  hnode.space.params.clear()
  node.parent = None
  for childnr, child in enumerate(children):
    assert isinstance(child, Node)
    node.attachChild(child, str(childnr), childnr)    
  return node
