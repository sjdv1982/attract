def import_link(link, is_module, glob = None):
  link = str(link)
  rest = link
  last = None
  lastmod = None
  pos = 0
  while 1:
    dot = link.find(".",pos+1)
    if dot == -1: 
      rest = link[pos:]
      break
    last = link[:dot]
    lastmod_old = lastmod
    lastmod = __import__(last)
    if lastmod is lastmod_old:
      lastmod = getattr(lastmod, link[pos:dot])
    pos = dot + 1
  if is_module:
    __import__(link)
  elif last is None:
    assert glob is not None
    return glob[rest]
  else:
    return getattr(lastmod, rest)

