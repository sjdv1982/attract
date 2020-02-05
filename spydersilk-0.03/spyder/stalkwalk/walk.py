def walk(node, cb_start, cb_end, cb_child_start, cb_child_end):
  cb_start(node)
  for childname in node.getChildren():
    child = node[childname]
    cb_child_start(childname)
    walk(child, cb_start, cb_end, cb_child_start, cb_child_end)
    cb_child_end(childname)
  cb_end(node)  