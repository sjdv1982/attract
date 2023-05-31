# Copyright 2007, Sjoerd de Vries
# This file is part of the Spyder module: "tarantula" 
# For licensing information, see LICENSE.txt 

import ctypes
from math import *

def calc_normal(vertices):
    normal = Coordinate(0,0,0)
    for nr in range(len(vertices)):
      lastnr = nr - 1
      if nr == 0: lastnr = len(vertices) - 1
      v1 = vertices[lastnr]
      v2 = vertices[nr]
      normal.x += v1.y * v2.z - v1.z * v2.y
      normal.y += v1.z * v2.x - v1.x * v2.z
      normal.z += v1.x * v2.y - v1.y * v2.x
    return normal.normalize()

def renderPolygon(a, *args, **kargs):
    vertices = []
    for v in a.vertices:
      #x = v.x * a.axis.x.x + v.y * a.axis.y.x + v.z * a.axis.z.x + a.axis.origin.x
      #y = v.x * a.axis.x.y + v.y * a.axis.y.y + v.z * a.axis.z.y + a.axis.origin.y
      #z = v.x * a.axis.x.z + v.y * a.axis.y.z + v.z * a.axis.z.z + a.axis.origin.z
      #vertices.append(Coordinate(x,y,z))
      vertices.append(Coordinate(v*a.axis))
      
    p = []
    for l in vertices:
        p += (l.x, l.y, l.z)
    pointlisttype = ctypes.c_double * len(p)
    pp = pointlisttype(*p)
    colortype = ctypes.c_ubyte * 3
    mat = get_material(a.material)
    color = colortype(mat.color.r, mat.color.g, mat.color.b)
    normal = calc_normal(vertices)
    normals = a.vertexnormals
    if normals == None:
      normals = []
      for nr in range(len(vertices)):
        normals += (normal.x, normal.y, normal.z)
    normals = pointlisttype(*normals)
    __render(len(vertices), pp, normals, color)

def prepare_renderObject(o):
    p = []
    for l in o.vertices:
      p += (l.x, l.y, l.z)
    pointlisttype = ctypes.c_double * len(p)
    points = pointlisttype(*p)
    points = ctypes.cast(points, ctypes.POINTER(ctypes.c_double))
    colors = []
    lengths = []
    indices = []
    for p in o.faces:
      if p.normal == None:
        verts = []
        for i in p.vertices: verts.append(o.vertices[i])

        p.normal = calc_normal(verts)
      mat = p.material
      if mat == None: mat = get_material(o.material)
      colors += (mat.color.r,mat.color.g,mat.color.b)
      lengths.append(len(p.vertices))
      indices += p.vertices
    colorstype = ctypes.c_ubyte * len(colors)
    colors = colorstype(*colors)
    colors = ctypes.cast(colors, ctypes.POINTER(ctypes.c_ubyte))
    lengthstype = ctypes.c_int * len(lengths)
    lengths = lengthstype(*lengths)
    lengths = ctypes.cast(lengths, ctypes.POINTER(ctypes.c_int))
    normalstype = ctypes.c_double * (len(indices) * 3)    
    indicestype = ctypes.c_int * len(indices)
    indices = indicestype(*indices)
    indices = ctypes.cast(indices, ctypes.POINTER(ctypes.c_int))
    axis = o.axis.matrix()
    axis = ctypes.cast(axis, ctypes.POINTER(ctypes.c_double))

    normals = []
    for p in o.faces:
      norms = p.vertexnormals
      if norms == None:
        norms = []
        if o.lighting == "flat":
          for pp in p.vertices:
            norms += [p.normal.x,p.normal.y,p.normal.z]
        elif o.lighting == "smooth":
          for pp in p.vertices:
            avgnormal = Coordinate(0,0,0)
            for ppp in o.faces:
              if pp in ppp.vertices:
                avgnormal += ppp.normal
            avgnormal = avgnormal.normalize()
            norms += [avgnormal.x,avgnormal.y,avgnormal.z]
      normals += norms
    normals = normalstype(*normals)
    normals = ctypes.cast(normals, ctypes.POINTER(ctypes.c_double))
    #print len(pp), len(colors), len(lengths), len(indices)
    return [ctypes.c_int(len(o.vertices)), points, ctypes.c_int(len(o.faces)), lengths, indices, colors, axis,normals]


def renderObject(o, show=True):
  #print [v for v in o.__conversionstack__ if type(v) == int]
  ptrs = prepare_renderObject(o)
  ptrs.append(show==False)
  displaylist = __renderObject(*ptrs)
  return (displaylist, o.axis)

def renderObjectArray(a, show=True):
  ret = []
  onr = 0
  if len(a) > 500:
    print("Please wait patiently while Spyder sends a list of %d objects..." % len(a))
  ptrs = []
  axes = []
  for o in a:
    o.__conversionstack__ = list(a.__conversionstack__)
    onr += 1
    if not onr % 500: print(onr, "/", len(a))
    ptrs.append(prepare_renderObject(o))
    axes.append(o.axis)
  newptrs = [len(ptrs)]    
  for nr in range(len(ptrs[0])):
    p = [a[nr] for a in ptrs]
    p = (type(p[0]) * len(p))(*p)
    newptrs.append(p)
  newptrs.append(show==False)
  ret0 = __renderObjectArray(*newptrs)
  ret = range(ret0,ret0+len(ptrs))
  return zip(ret,axes)

def make_displaylist(l):
  lists = []
  axes0 = []
  counter = 1
  while 1:
    newl = []
    counter2 = 0
    for ll in l:
      counter2 +=1
      #print counter, counter2
      #print ll
      if isinstance(ll[0],int):
        if int(ll[0]) > -1:
          lists.append(int(ll[0]))
          axes0.append(ll[1])      
      else: newl += ll
    if len(newl) == 0: break
    l = newl
  nrlists = len(lists)
  liststype = ctypes.c_int * nrlists
  lists = liststype(*lists)
  axes = []
  for a in axes0:
    axes += [a.x.x,a.x.y,a.x.z,0,
            a.y.x,a.y.y,a.y.z,0,
            a.z.x,a.z.y,a.z.z,0,
            a.origin.x, a.origin.y, a.origin.z, 1]  
  axestype = ctypes.c_double*(16*nrlists)
  axes = axestype(*axes)
  displaylist = __assembleDisplayList(nrlists,lists,axes)
  return displaylist
  
    
def axissystem_to_matrix(a):
  import ctypes
  axistype = ctypes.c_double * 16
  axis = [a.x.x,a.x.y,a.x.z,0,
        a.y.x,a.y.y,a.y.z,0,
        a.z.x,a.z.y,a.z.z,0,
        a.origin.x, a.origin.y, a.origin.z, 1]
  return ctypes.cast(axistype(*axis), ctypes.POINTER(ctypes.c_double))
    

def renderDisplayList(d,show=True):
  import ctypes
  axis = d.axis.matrix()
  axestype =  ctypes.POINTER(ctypes.c_double) * 1
  axes = axestype(axis)
  __renderList(d.displaylist, axes,1,show==True)
  return (d.displaylist, d.axis)

def renderMultiDisplayList(md,show=True):
  import ctypes
  if len(md.instances) == 0: return None
  axestype =  ctypes.POINTER(ctypes.c_double) * len(md.instances)  
  axes = []
  for a in md.instances: axes.append(a.matrix())
  axes = axestype(*axes)
  ret = __renderList(md.displaylist, axes,len(md.instances),show==True)
  return (ret, AxisSystem())

    #print list(lengths)
def makeBlock(b):
  if b.pivot == "center":
    minmax = (-1,1)
    sc = Coordinate(b.dimensions / 2)
  elif b.pivot == "corner":
    minmax = (0,1)
    sc = Coordinate(b.dimensions)
  vertices = []
  for n in minmax:
      for nn in minmax:
          for nnn in minmax:
              vertices.append(Coordinate(sc.x*n, sc.y*nn,sc.z*nnn))
  vertices = CoordinateArray(vertices)
  faces = []
  #indices = ((0,1,3,2), (4,5,7,6), (0,4,5,1), (1,5,7,3), (3,7,6,2), (2,6,4,0))
  indices = ((0,1,3,2), (4,6,7,5), (0,4,5,1), (1,5,7,3), (3,7,6,2), (2,6,4,0))
  for n in range(0,6):
    vertices0 = indices[n]
    p = Face3D(vertices=vertices0)
    faces.append(p)
  o = Object3D(vertices=vertices, faces=faces,material=b.material, axis=b.axis)
  return o


def drawCircle(c, segments):
    r = c.radius
    vertices = []
    a = Vector(1,0,0)
    if a * c.normal > .999: a = Vector(0,1,0)
    inplane1 = a ^ c.normal
    inplane1 = inplane1.normalize()    
    inplane2 = inplane1 ^ c.normal
    for n in range(0,segments):
      ang = float(n)/segments * 2 * pi
      cosang = cos(ang) * r
      sinang = sin(ang) * r
      x = cosang * inplane1.x + sinang * inplane2.x
      y = cosang * inplane1.y + sinang * inplane2.y
      z = cosang * inplane1.z + sinang * inplane2.z
      vertices.append(Coordinate(x,y,z) + c.origin)
    p = Face3D(material=c.material,vertices=range(0,segments))
    o = Object3D(vertices=vertices, faces=(p,),material=c.material)
    return o

def drawCylinder(c, segments):
    r = c.radius
    z = c.axis.z.normalize()
    vertices = []
    a = Vector(0,1,0)
    if z.y > .999: a = Vector(1,0,0)
    inplane1 = a ^ z
    inplane2 = z  ^ inplane1
    for n in range(0,segments):
      ang = float(n)/segments * 2 * pi
      cosang = cos(ang) * r
      sinang = sin(ang) * r
      vertices.append(Coordinate(cosang,sinang,-0.5*c.height))
      vertices.append(Coordinate(cosang,sinang,0.5*c.height))      
    faces = []
    v = range(0,2*segments,2)    
    v.reverse()
    p = Face3D(vertices=v)
    faces.append(p)
    v = range(1,2*segments,2)
    p = Face3D(vertices=v)
    faces.append(p)        
    for n in range(0,segments):
        nextn = n + 1
        if nextn == segments: nextn = 0
        p = Face3D(vertices=(2*n+1, 2*n, 2*nextn, 2*nextn+1))            
        faces.append(p)
    axis=AxisSystem(x=inplane1,y=inplane2,z=c.axis.z,origin=c.axis.origin)
    axis.x = axis.x * c.axis.x.size()
    axis.y = axis.y * c.axis.y.size()
    o = Object3D(vertices=vertices, faces=faces,axis=axis, material=c.material)
    return o



def drawCircle16(c): return drawCircle(c,16)
def drawCircle32(c): return drawCircle(c,32)
def drawCircle64(c): return drawCircle(c,64)

def drawCylinder16(c): return drawCylinder(c,16)
def drawCylinder32(c): return drawCylinder(c,32)
def drawCylinder64(c): return drawCylinder(c,64)


def untieObject(o):
    polygons = []
    for p in o.faces:
        vertices = []
        for v in p.vertices:
            vertices.append(o.vertices[v])
        mat = p.material
        if mat == None: mat = o.material
        pp = Polygon(vertices, mat, o.axis, p.texturecoords, p.normal, p.vertexnormals)
        polygons.append(pp)
    return PolygonArray(polygons)
    
def renderMultiInstance(mi, show=True):
  dlist = make_displaylist(mi.object.show(show=False))
  axes = []
  import ctypes
  for a in mi.instances:
    p = a.matrix()
    axes.append(p)
  if len(axes) > 0:
    axestype =  ctypes.POINTER(ctypes.c_double) * len(axes)
    axes2 = axestype(*axes)
    ret = __renderList(dlist, axes2,len(axes),int(show==True))
  else:
    ret = -1
  return ret, AxisSystem()
  
_materials = {}

def renderObjectGroup(og, show=True):
  dlist = make_displaylist(og.group.show(show=False))
  import ctypes
  axestype =  ctypes.POINTER(ctypes.c_double) * 1
  axes = axestype(og.axis.matrix())  
  ret = __renderList(dlist, axes,1,int(show==True))
  return ret, AxisSystem()
  
def generate_material_name():
  n = 0
  while 1:
      n += 1
      s = "Material" + str(n)
      if s not in _materials: break
  return s

def register_newmaterial(mat,*args,**kargs):
  assert mat.typename() == "NewMaterial"
  if mat.name == None: mat.name = generate_material_name()
  _materials[mat.name] = mat
  return mat.name

def get_material(name):
  if name not in _materials: raise KeyError("material %s has not been defined" % name) 
  return _materials[name]
  
