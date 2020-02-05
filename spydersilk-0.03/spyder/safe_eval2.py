#----------------------------------------------------------------------
# I, Babar K. Zafar, the author or of this code dedicate any and all
# copyright interest in this code to the public domain. I make this
# dedication for the benefit of the public at large and to the
# detriment of our heirs and successors. I intend this dedication to
# be an overt act of relinquishment in perpetuity of all present and
# future rights this code under copyright law.
#
# Version 0.1 / May 27 2006
# Jan - Feb 2007: adjustments by SJdV 
# - implementing a suggestion by Giovanni Bajo
#   to enable "import safe_eval": (__builtins__ => __builtin__)
# - re-enabling exceptions
# - disabling context checking
# - disabling execution itself, only checking for unsafe code
#----------------------------------------------------------------------

import __builtin__
import inspect, compiler.ast
#import thread, time

#----------------------------------------------------------------------
# Module globals.
#----------------------------------------------------------------------

# Toggle module level debugging mode.
DEBUG = False

# List of all AST node classes in compiler/ast.py.
all_ast_nodes = \
    [name for (name, obj) in inspect.getmembers(compiler.ast)
     if inspect.isclass(obj) and issubclass(obj, compiler.ast.Node)]

# List of all builtin functions and types (ignoring exception classes).
all_builtins = \
    [name for (name, obj) in inspect.getmembers(__builtin__)
     if inspect.isbuiltin(obj) or (inspect.isclass(obj) and \
                                   not issubclass(obj, Exception))]

#----------------------------------------------------------------------
# Utilties.
#----------------------------------------------------------------------

def classname(obj):
    return obj.__class__.__name__

def is_valid_ast_node(name):
    return name in all_ast_nodes

def is_valid_builtin(name):
    return name in all_builtins

def get_node_lineno(node):
    return (node.lineno) and node.lineno or 0
       
#----------------------------------------------------------------------
# Restricted AST nodes & builtins.
#----------------------------------------------------------------------

# Deny evaluation of code if the AST contain any of the following nodes:
unallowed_ast_nodes = [
#   'Add', 'And',
#   'AssAttr', 'AssList', 'AssName', 'AssTuple',
#   'Assert', 'Assign', 'AugAssign',
    'Backquote',
#   'Bitand', 'Bitor', 'Bitxor', 'Break',
#   'CallFunc', 'Class', 'Compare', 'Const', 'Continue',
#   'Decorators', 'Dict', 'Discard', 'Div',
#   'Ellipsis', 'EmptyNode',
    'Exec',
#   'Expression', 'FloorDiv',
#   'For',
    'From',
#   'Function',
#   'GenExpr', 'GenExprFor', 'GenExprIf', 'GenExprInner',
#   'Getattr', 'Global', 'If',
    'Import',
#   'Invert',
#   'Keyword', 'Lambda', 'LeftShift',
#   'List', 'ListComp', 'ListCompFor', 'ListCompIf', 'Mod',
#   'Module',
#   'Mul', 'Name', 'Node', 'Not', 'Or', 'Pass', 'Power',
#   'Print', 'Printnl',
#    'Raise',
#    'Return', 'RightShift', 'Slice', 'Sliceobj',
#   'Stmt', 'Sub', 'Subscript',
#   'TryExcept', 'TryFinally',
#   'Tuple', 'UnaryAdd', 'UnarySub',
#   'While','Yield'
]

# Deny evaluation of code if it tries to access any of the following builtins:
unallowed_builtins = [
    '__import__',
#   'abs', 'apply', 'basestring', 'bool', 'buffer',
#   'callable', 'chr', 'classmethod', 'cmp', 'coerce',
    'compile',
#   'complex',
    'delattr',
#   'dict',
    'dir',
#   'divmod', 'enumerate',
    'eval', 'execfile', 'file',
#   'filter', 'float', 'frozenset',
    'globals', 
    'getattr', 'hasattr', 
#    'hash', 'hex', 'id',
    'input',
#   'int', 'intern', 'isinstance', 'issubclass', 'iter',
#   'len', 'list',
    'locals',
#   'long', 'map', 'max', 'min', 'object', 'oct',
    'open',
#   'ord', 'pow', 'property', 'range',
    'raw_input',
#   'reduce',
    'reload',
#   'repr', 'reversed', 'round', 'set',
    'setattr',
#   'slice', 'sorted', 'staticmethod',  'str', 'sum', 'super',
#   'tuple', 'type', 'unichr', 'unicode',
    'vars',
#    'xrange', 'zip'
]

for ast_name in unallowed_ast_nodes:
    assert(is_valid_ast_node(ast_name))
for name in unallowed_builtins:
    if name in ('input', 'raw_input'): continue #not valid builtins in IPython
    assert(is_valid_builtin(name)), name 

def is_unallowed_ast_node(kind):
    return kind in unallowed_ast_nodes

def is_unallowed_builtin(name):
    return name in unallowed_builtins

#----------------------------------------------------------------------
# Restricted attributes.
#----------------------------------------------------------------------

# In addition to these we deny access to all lowlevel attrs (__xxx__).
unallowed_attr = [
    'im_class', 'im_func', 'im_self',
    'func_code', 'func_defaults', 'func_globals', 'func_name',
    'tb_frame', 'tb_next',
    'f_back', 'f_builtins', 'f_code', 'f_exc_traceback',
    'f_exc_type', 'f_exc_value', 'f_globals', 'f_locals']

    
def is_unallowed_attr(name):
    if (name[:2] == '__'):
         return True
    else: return (name in unallowed_attr)

#----------------------------------------------------------------------
# SafeEvalVisitor.
#----------------------------------------------------------------------

class SafeEvalError(object):
    """
    Base class for all which occur while walking the AST.

    Attributes:
      errmsg = short decription about the nature of the error
      lineno = line offset to where error occured in source code
    """
    def __init__(self, errmsg, lineno):
        self.errmsg, self.lineno = errmsg, lineno
    def __str__(self):
        return "line %d : %s" % (self.lineno, self.errmsg)

class SafeEvalASTNodeError(SafeEvalError):
    "Expression/statement in AST evaluates to a restricted AST node type."
    pass
class SafeEvalBuiltinError(SafeEvalError):
    "Expression/statement in tried to access a restricted builtin."
    pass
class SafeEvalAttrError(SafeEvalError):
    "Expression/statement in tried to access a restricted attribute."
    pass

class SafeEvalVisitor(object):
    """
    Data-driven visitor which walks the AST for some code and makes
    sure it doesn't contain any expression/statements which are
    declared as restricted in 'unallowed_ast_nodes'. We'll also make
    sure that there aren't any attempts to access/lookup restricted
    builtin declared in 'unallowed_builtins'. By default we also won't
    allow access to lowlevel stuff which can be used to dynamically
    access non-local envrioments.

    Interface:
      walk(ast) = validate AST and return True if AST is 'safe'

    Attributes:
      errors = list of SafeEvalError if walk() returned False

    Implementation:
    
    The visitor will automatically generate methods for all of the
    available AST node types and redirect them to self.ok or self.fail
    reflecting the configuration in 'unallowed_ast_nodes'. While
    walking the AST we simply forward the validating step to each of
    node callbacks which take care of reporting errors.
    """

    def __init__(self):
        "Initialize visitor by generating callbacks for all AST node types."
        self.errors = []
        for ast_name in all_ast_nodes:
            # Don't reset any overridden callbacks.
            if getattr(self, 'visit' + ast_name, None): continue
            if is_unallowed_ast_node(ast_name):
                setattr(self, 'visit' + ast_name, self.fail)
            else:
                setattr(self, 'visit' + ast_name, self.ok)

    def walk(self, ast):
        "Validate each node in AST and return True if AST is 'safe'."
        self.visit(ast)
        return self.errors == []
        
    def visit(self, node, *args):
        "Recursively validate node and all of its children."
        fn = getattr(self, 'visit' + classname(node))
        if DEBUG: self.trace(node)
        fn(node, *args)
        for child in node.getChildNodes():
            self.visit(child, *args)

    def visitName(self, node, *args):
        "Disallow any attempts to access a restricted builtin/attr."
        name = node.getChildren()[0]
        lineno = get_node_lineno(node)
        if is_unallowed_builtin(name):
            self.errors.append(SafeEvalBuiltinError( \
                "access to builtin '%s' is denied" % name, lineno))
        elif is_unallowed_attr(name):
            self.errors.append(SafeEvalAttrError( \
                "access to attribute '%s' is denied" % name, lineno))
    
    def visitAssAttr(self, node, *args):
      return self.visitGetattr(node, *args)

    def visitGetattr(self, node, *args):
        "Disallow any attempts to access a restricted attribute."
        name = node.attrname
        lineno = get_node_lineno(node)
        if is_unallowed_attr(name):
            self.errors.append(SafeEvalAttrError( \
                "access to attribute '%s' is denied" % name, lineno))
            
    def ok(self, node, *args):
        "Default callback for 'harmless' AST nodes."
        pass
    
    def fail(self, node, *args):
        "Default callback for unallowed AST nodes."
        lineno = get_node_lineno(node)
        self.errors.append(SafeEvalASTNodeError( \
            "execution of '%s' statements is denied" % classname(node),
            lineno))

    def trace(self, node):
        "Debugging utility for tracing the validation of AST nodes."
        print classname(node)
        for attr in dir(node):
            if attr[:2] != '__':
                print ' ' * 4, "%-15.15s" % attr, getattr(node, attr)

#----------------------------------------------------------------------
# Safe 'eval' replacement.
#----------------------------------------------------------------------

class SafeEvalException(Exception):
    "Base class for all safe-eval related errors."
    pass

class SafeEvalCodeException(SafeEvalException):
    """
    Exception class for reporting all errors which occured while
    validating AST for source code in safe_eval().

    Attributes:
      code   = raw source code which failed to validate
      errors = list of SafeEvalError
    """
    def __init__(self, code, errors):
        self.code, self.errors = code, errors
    def __str__(self):
        return '\n'.join([str(err) for err in self.errors])

        
def safe_eval(code):
    """
    Validate source code and make sure it contains no unauthorized
    expression/statements as configured via 'unallowed_ast_nodes' and
    'unallowed_builtins'. By default this means that code is not
    allowed import modules or access dangerous builtins like 'open' or
    'eval'. The following exception will be raised on errors:

      if code is didn't validate and is considered 'unsafe' = 
        SafeEvalCodeException

    """

    try:
      ast = compiler.parse(code)
    except SyntaxError:
      ast = None
      pass #will crash anyway later, and we will know the line number
    except:
      print code
      raise
    if ast != None:
      checker = SafeEvalVisitor()

      if checker.walk(ast) == False:
        raise SafeEvalCodeException(code, checker.errors)

"""
if __name__ == "__main__":
    s = "a = Integer(3)\n"
    print safe_eval(s)
"""    


