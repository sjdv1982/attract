from __future__ import print_function, absolute_import

names = ["QtGui", "QtCore", "QtDeclarative", "QtWebKit"]

import sys
try:
  if "--pyqt" in sys.argv: raise ImportError
  import PySide
  _qt = "PySide"
  names.append("QtUiTools")
except ImportError:
  try:
    import PyQt4
    import sip
    for name in ["QDate", "QDateTime", "QString", "QTextStream", "QTime", "QUrl", "QVariant"]:
      sip.setapi(name, 2)    
    _qt = "PyQt4"
  except ImportError:
    raise ImportError("You must have either PySide OR PyQt4 installed")


_qtmod = __import__(_qt, globals(), locals(), names)
for name in names:
  mod = getattr(_qtmod, name)
  globals()[name] = mod
  sys.modules[__name__ + "." + name] = mod

if _qt == "PyQt4":
  import PyQt4.uic
  import imp
  class QUiLoader(object):
    @staticmethod
    def load(arg):
      if isinstance(arg,QtCore.QBuffer):
        buf = arg
        if sys.version_info[0] == 3:
          s = bytes(buf.data())
          s = s.decode() +"\n"
        else:
          s = str(buf.data())+"\n"
          s = unicode(s)
        try:        
          from io import StringIO
        except ImportError:
          from cStringIO import StringIO        
        strio = StringIO(s)
        return PyQt4.uic.loadUi(strio)
      elif isinstance(arg,QtCore.QFile):
        arg.open(QtCore.QFile.ReadOnly)
        return PyQt4.uic.loadUi(arg)
      else: 
        raise TypeError(type(arg))
  QtUiTools = imp.new_module(__name__ + ".QtUiTools")
  QtUiTools.QUiLoader = QUiLoader
  sys.modules[__name__ + ".QtUiTools"] = QtUiTools
elif _qt == "PySide":
  QtCore.pyqtSignal = QtCore.Signal
  QtCore.pyqtSlot = QtCore.Slot  

  
  
