class parameters(object): 
  def __init__(self, **args):
    for a in args: setattr(self, a, args[a])

def reserve_slot(slot):
  import bee
  from bee.types import stringtupleparser
  from dragonfly.canvas import box2d
  c = bee.configure(slot.canvasname)
  b = slot.box
  box = box2d(b.x, b.y, b.sizex, b.sizey, b.mode)
  params = None
  col = slot.color
  if col is not None:
    params = parameters(color=(col.r,col.g,col.b,col.a))
  c.reserve (
   slot.slotname,
   stringtupleparser(slot.slottype),
   box=box,
   parameters=params,
  )
  return c
