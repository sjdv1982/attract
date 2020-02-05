from .xml import _group_generators, wrappings_regen_min, wrappings_regen_max, update
def resolve_elegroup(
 elegroup,
 elegroups,
 xml,
 t,
 curr_wrappings, 
 done_membernames,
 mxmls,
 mxmldata,
 mwrappings,
 curr_membernames,
 curr_members,
 mcomplexities,
 groupindex,
 name,
 n,
 objs,
 forms,
 wrappings,
 complexity,
 gen_group,
):
  g = find_group(elegroups, elegroup.id)
  if g is not None and g.id not in done_membernames:
    #resolve the parent group instead
    return resolve_elegroup(
            g,
            elegroups,
            xml,
            t,
            curr_wrappings, 
            done_membernames,
            mxmls,
            mxmldata,
            mwrappings,
            curr_membernames,
            curr_members,
            mcomplexities,
            groupindex,
            name,
            n,
            objs,
            forms,
            wrappings,
            complexity,
            gen_group,
           )

  done_membernames.add(elegroup.id)  
  elegroupids = [e.id for e in elegroups]
  gen = _group_generators[elegroup.type, complexity]    
  start = gen(name, n, objs, forms, wrappings, "start", elegroup)
  if start is not None: 
    update(xml, start, t)
    t.incr()
    curr_wrappings += 1
  for mem in elegroup.members:
    if mem not in curr_membernames and mem not in elegroupids:
      raise Exception("Unknown member %s in group element %s" % (mem, elegroup.id))
    pre = gen(name, n, objs, forms, wrappings, "pre", elegroup)
    if pre is not None: 
      update(xml, pre, t)
      t.incr()
      curr_wrappings += 1
    while 1:
      if mem in elegroupids:

        g = [e for e in elegroups if e.id == mem][0]
        curr_wrappings, groupindex = resolve_elegroup(
         g,
         elegroups,
         xml,
         t,
         curr_wrappings, 
         done_membernames,
         mxmls,
         mxmldata,
         mwrappings,
         curr_membernames,
         curr_members,
         mcomplexities,
         groupindex,
         name,
         n,
         objs,
         forms,
         wrappings,
         complexity,
         gen_group,
        )      
        break

      for mxml, mdat, membername, member, mcomp in \
       zip(mxmls, mxmldata, curr_membernames, curr_members, mcomplexities):
        if membername != mem: continue
        curr_wrappings, groupindex = resolve_member(
         membername,
         elegroups,
         xml,
         t,
         curr_wrappings, 
         done_membernames,
         mxmls,
         mxmldata,
         mwrappings,
         curr_membernames,
         curr_members,
         mcomplexities,
         groupindex,
         name,
         n,
         objs,
         forms,
         wrappings,
         complexity,
         gen_group,
        )
      post = gen(name, n, objs, forms, wrappings, "post", elegroup)
      if post is not None: 
        update(xml, post, t)
        t.decr()
        curr_wrappings -= 1

      break  
  end = gen(name, n, objs, forms, wrappings, "end", elegroup)
  if end is not None: 
    update(xml, end, t)
    t.decr()
    curr_wrappings -= 1

  return curr_wrappings, groupindex    
      

def find_group(elegroups, membername):
  ret = None
  for g in elegroups:
    if membername in g.members: 
      if ret is not None:
        raise Exception("membername %s is a member of multiple group elements" % membername)
      ret = g
  return ret     
  

def resolve_member(
 membername,
 elegroups,
 xml,
 t,
 curr_wrappings, 
 done_membernames,
 mxmls,
 mxmldata,
 mwrappings,
 curr_membernames,
 curr_members,
 mcomplexities,
 groupindex,
 name,
 n,
 objs,
 forms,
 wrappings,
 complexity,
 gen_group,
):
  from .xml import _write_xml
  for mxml, mdat, membername0, member, mcomp in \
   zip(mxmls, mxmldata, curr_membernames, curr_members, mcomplexities):
    if membername0 != membername: continue
    if membername in done_membernames: continue
    g = find_group(elegroups, membername)
    if g is not None and g.id not in done_membernames:
      curr_wrappings, groupindex = resolve_elegroup(
       g,
       elegroups,
       xml,
       t,
       curr_wrappings, 
       done_membernames,
       mxmls,
       mxmldata,
       mwrappings,
       curr_membernames,
       curr_members,
       mcomplexities,
       groupindex,
       name,
       n,
       objs,
       forms,
       wrappings,
       complexity,
       gen_group,
      )      
      continue
    done_membernames.add(membername)
    if mcomp > 1 and curr_wrappings != mwrappings:
      if max(mwrappings, curr_wrappings) >= wrappings_regen_min \
        and min(mwrappings, curr_wrappings) <= wrappings_regen_max:
         #wrapping level has changed, regenerate XML...
         mobjs, mforms, newprename = mdat
         mxml, dmmy = _write_xml(mobjs, mforms, membername,newprename,curr_wrappings)
    for l in mxml[0]:
      xml[0].append(t.tabstr + l)
    for l in mxml[1]:
      xml[1].append(t.tabstr2 + l)
  
  return curr_wrappings, groupindex    
  
def resolve_group(
 group,
 elegroups,
 xml,
 t,
 curr_wrappings, 
 done_membernames,
 mxmls,
 mxmldata,
 mwrappings,
 curr_membernames,
 curr_members,
 mcomplexities,
 groupindex,
 name,
 n,
 objs,
 forms,
 wrappings,
 complexity,
 gen_group,
):
  groupindex += 1
  start, end = None, None
  if gen_group is not None:
    start = gen_group(name, n, objs, forms, wrappings, group, groupindex, "start")
  if start is not None: 
    update(xml, start, t)
    t.incr()
    curr_wrappings += 1
  for membername2, member2 in zip(curr_membernames, curr_members):

    if membername2 in done_membernames: continue    
    if hasattr(member2, "group") and member2.group == group:
      pre = None
      if gen_group is not None:
        pre = gen_group(name, n, objs, forms, wrappings, group, groupindex, "pre")
      if pre is not None: 
        update(xml, pre, t)
        t.incr()
        curr_wrappings += 1

      curr_wrappings, groupindex = resolve_member (
        membername2,
        elegroups,
        xml,
        t,
        curr_wrappings, 
        done_membernames,
        mxmls,
        mxmldata,
        mwrappings,
        curr_membernames,
        curr_members,
        mcomplexities,
        groupindex,
        name,
        n,
        objs,
        forms,
        wrappings,
        complexity,
        gen_group,
      )
      post = None
      if gen_group is not None:
        post = gen_group(name, n, objs, forms, wrappings, group, groupindex, "post")
      if post is not None: 
        update(xml, post, t)
        t.decr()
        curr_wrappings -= 1

  if gen_group is not None:
    end = gen_group(name, n, objs, forms, wrappings, group, groupindex, "end")
  if end is not None: 
    update(xml, end, t)
    t.decr()
    curr_wrappings -= 1
  
  return curr_wrappings, groupindex

def resolve(
 xml,
 t,
 curr_wrappings, 
 done_membernames,
 mxmls,
 mxmldata,
 mwrappings,
 curr_membernames,
 curr_members,
 mcomplexities,
 groupindex,
 name,
 n,
 objs,
 forms,
 wrappings,
 complexity,
 gen_group,
 toplevel = False
):
  form = forms[0]
  elegroups = form._groups
  elegroupmembers = set()
  for e in elegroups:
    for m in e.members: elegroupmembers.add(m)
  for membername, member in zip(curr_membernames, curr_members):
    if hasattr(member, "group") and member.group is not None:
      if membername in elegroupmembers:
        raise Exception("member %s has a .group attribute and is also a member of a group element" % membername)

  for mxml, membername, member in zip(mxmls, curr_membernames, curr_members):
    if membername in done_membernames: continue    
    if hasattr(member, "group") and member.group is not None:
      group = member.group
      curr_wrappings, groupindex = resolve_group (
       group,
       elegroups,
       xml,
       t,
       curr_wrappings, 
       done_membernames,
       mxmls,
       mxmldata,
       mwrappings,
       curr_membernames,
       curr_members,
       mcomplexities,
       groupindex,
       name,
       n,
       objs,
       forms,
       wrappings,
       complexity,
       gen_group,
      )
    else:      
      curr_wrappings, groupindex = resolve_member (
        membername,
        elegroups,
        xml,
        t,
        curr_wrappings, 
        done_membernames,
        mxmls,
        mxmldata,
        mwrappings,
        curr_membernames,
        curr_members,
        mcomplexities,
        groupindex,
        name,
        n,
        objs,
        forms,
        wrappings,
        complexity,
        gen_group,
      )

  #TODO: also for groups, elegroups!!
  newxml = []
  memberpos = 0
  for l in xml[0]:
    if toplevel == False:
      if l.find("!!NEW-MEMBERPOS!!") > -1:
        memberpos += 1
        l = l.replace("!!NEW-MEMBERPOS!!", str(memberpos-1))
      elif l.find("!!MEMBERPOS!!") > -1:
        l = l.replace("!!MEMBERPOS!!", str(memberpos-1))
    newxml.append(l)

  xml[0][:] = newxml
