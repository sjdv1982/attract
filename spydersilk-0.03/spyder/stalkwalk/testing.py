from __future__ import print_function

def print_grouping(root):
  name = []
  def cb_start(node):
    if node.name is None: return
    name.append(node.name)   
    space = node.space
    if space.topgroup is not None:
      fullname = "-".join([str(n) for n in name])
      print("GROUPING", fullname)
      def print_group(group, groupname):
        if groupname is not None: print("GROUP", groupname)
        for member, membertype in group.members:
          if membertype == "group":
            print_group(space.groups[member], member)
          else:
            print(member)
        if groupname is not None: print("/GROUP", groupname)
      print_group(space.topgroup, None)
      print("/GROUPING", fullname)
      
  def cb_end(node):
    if len(name): name[:] = name[:-1]
    
  def cb_none(*args, **kwargs): pass  
  from .walk import walk
  walk(root, cb_start, cb_end, cb_none, cb_none)
 
def print_tree(root): 
  indent = []
  def cb_start(node):
    if len(indent) <= 999: 
      print(len(indent)*"  " + "START", node.space.wrapping, node.space.complexity, node.name, node.nodetype, node.space.formpath)
    indent.append(0)
  def cb_end(node):
    indent.pop(0)
    #if len(indent) <= 2: 
    #  print(len(indent)*"  " + "END", node.name, node.nodetype)
  def cb_none(*args, **kwargs): pass    
  from .walk import walk
  walk(root, cb_start, cb_end, cb_none, cb_none)
