from . import Handler

class Group(object):
  def __init__(self, name, type, elegroup = None):
    self.name = name
    self.type = type
    if self.type == "elegroup":
      assert elegroup is not None
      self.elegroup = elegroup
    else:
      assert elegroup is None
    self.members = []
  def add_member(self, member, membertype):
    assert isinstance(member, str) or isinstance(member, int), member
    assert membertype in ("member", "submember", "group"), membertype
    self.members.append((member, membertype))

def determine_groups(message, node):
  """
  Make an initial organization of all members into groups
   based on the .group and ._groups (elegroups) attributes
  After this, the groups must be pruned:
  - The topgroup using grouped_members and removed_members
  - All other groups using removed_members
  """
  if not node.getChildren(): return Handler.MATCH
  space = node.space
  
  #Create four collections
  #topgroup: initially the memberlist, later also with groups inserted
  #groups: a flat dict of groups
  #groupnames: the names of the groups in order
  #grouped_members: members that are now inside a group
  #removed_members: members removed by the parent to a higher-level group
    
  curr_membernames = space.curr_membernames
  assert curr_membernames
  membernames = []
  if space.memberorder is not None:
    for mnr, m in enumerate(space.memberorder):
      assert m not in space.memberorder[mnr+1:], m
      if m not in curr_membernames: continue
      membernames.append(m)  
  for membername in curr_membernames:
    if membername not in membernames:
      membernames.append(membername)
  
  space.topgroup = Group(None, None)
  for membername in membernames:
    space.topgroup.add_member(membername, "member")
  space.groups = {}
  space.groupnames = []
  space.grouped_members = []
  if space.removed_members is None: space.removed_members = []
  
  offset = 0
  for nr,childname in enumerate(membernames):
    child = node[childname]
    if child.space.group is not None:
      groupname = child.space.group
      if groupname not in space.groups:
        space.groups[groupname] = Group(groupname, "group")
        space.groupnames.append(groupname)
        space.topgroup.members.insert(nr+offset, (groupname, "group"))
        offset += 1
      space.groups[groupname].add_member(childname, "member")  
      space.grouped_members.append(childname)      
    
  elegroupnames = [e.id for e in space.elegroups]  
  for e in space.elegroups:
    groupname = e.id
    assert groupname not in space.groups, groupname
    group = Group(groupname, "elegroup", e)
    space.groups[groupname] = group
    space.groupnames.append(groupname)
    for member in e.members:
      if member in membernames:
        mtype = "member"
      elif member.find(".") > -1 or member.find("[") > -1:        
        mtype = "submember"
      elif member in space.groups or member in elegroupnames:
        mtype = "group"
      else:
        raise ValueError("Unknown member", member, node.nodetype, node.name)
      group.add_member(member, mtype)
      if mtype == "member":
        assert member not in space.grouped_members, member
        space.grouped_members.append(member)              
      elif mtype == "submember":
        from .Node import decompose_attribute
        attr = decompose_attribute(member)
        node.sendMessage("grouping-remove-member", (attr,))
  
  for e in reversed(space.elegroups):
    groupname = e.id
    for ee in space.elegroups:
      if e is ee: continue
      if groupname in ee.members: break
    else:
      pos = 0
      if hasattr(e, "insert_at_member") and e.insert_at_member:
        for member in e.members:
          if member not in membernames: continue        
          for ind, t in enumerate(space.topgroup.members):
            tmember,tmtype = t
            if tmember == member and tmtype == "member":
              pos = ind
          else:
            continue
          break
      space.topgroup.members.insert(pos, (groupname, "group"))
    
  Handler.passdown(message, node)
  return Handler.MATCH

def prune_groups(message, node):
  if not node.getChildren(): return Handler.MATCH
  space = node.space
  for item in list(space.topgroup.members):
    if item[1] == "member" and \
     item[0] in space.removed_members or item[0] in space.grouped_members:
      space.topgroup.members.remove(item)
  for g in space.groups:
    group = space.groups[g]
    for item in list(group.members):
      if item[1] == "member" and item[0] in space.removed_members:
        group.members.remove(item)
  Handler.passdown(message, node)
  return Handler.MATCH
        
def make_groups(groupnodetype, message, node):
  if not node.getChildren(): return Handler.MATCH
  space = node.space
  topgroupmembers = set([v[0] for v in space.topgroup.members])
  groupindex = 0
  for pos, item in enumerate(space.topgroup.members):
    member, membertype = item
    if membertype != "group": continue
    if not space.groups[member].members: continue
    space.params.groups = space.groups
    space.params.groupindex = groupindex
    groupindex += 1
    space.params.formnode = node
    
    #determine real position: grouped-away members are still there as children...
    cpos2 = 0
    for cpos, childname in enumerate(list(node.getChildren())):
      if cpos2 == pos: break
      if childname in topgroupmembers: cpos2 += 1
    node.insertChild(groupnodetype, member, cpos)
    space.params.clear()
  Handler.passdown(message, node)
  return Handler.MATCH
  
def remove_member(message, node, attribute):
  assert isinstance(attribute, tuple), attribute
  assert len(attribute)
  for a in attribute:
    assert isinstance(a, str) or isinstance(a, int), a
  space = node.space  
  childname = attribute[0]
  if len(attribute) == 1:    
    if space.removed_members is None: space.removed_members = []
    assert childname in space.curr_membernames, (childname, space.curr_membernames)
    assert childname not in space.removed_members, childname
    space.removed_members.append(childname)
  else:
    node[childname].sendMessage(message, args = (attribute[1:],) )
  return Handler.MATCH

def init_group(node):
  params = node.parent.space.params
  
  space = node.space
  space.groupindex = params.groupindex
  assert isinstance(space.groupindex, int)
  space.groups = params.groups
  space.formnode = params.formnode
  
  space.group = space.groups[node.name]  
  if space.group.type == "elegroup":
    space.elegroup = space.group.elegroup
  groupindex = 0
  for member, membertype in space.group.members:
    if membertype == "group":
      space.params.groups = space.groups
      space.params.groupindex = groupindex
      groupindex += 1
      space.params.formnode = space.formnode
      node.addChild(node.nodetype, member)
      space.params.clear()
    else:
      node.addChild("dummy", member)
  return Handler.MATCH
  
def complete_group(message, node):  
  space = node.space
  group = space.group
  formnode = space.formnode
  for member, membertype in group.members:
    if membertype == "member":
      formnode.sendMessage("reparent", args = ((member,), node, member, None, None, True)) #keep_link
    elif membertype == "submember":
      from .Node import decompose_attribute
      attr = decompose_attribute(member)
      formnode.sendMessage("reparent", args = (attr, node, member, None, None, True)) #keep_link
  Handler.passdown(message, node)
  return Handler.MATCH
      