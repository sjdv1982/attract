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
#Jul 2009, adjustments by SJdV: making compatible with Python 3.0
#----------------------------------------------------------------------


import builtins
import inspect, ast
#import thread, time

#----------------------------------------------------------------------
# Module globals.
#----------------------------------------------------------------------

# Toggle module level debugging mode.
DEBUG = False

# List of all AST node classes in ast
import sys
all_ast_nodes = \
    [name for (name, obj) in inspect.getmembers(ast)
     if inspect.isclass(obj) and issubclass(obj, ast.AST)]

# List of all builtin functions and types (ignoring exception classes).
all_builtins = \
    [name for (name, obj) in inspect.getmembers(builtins)
     if inspect.isbuiltin(obj) or (inspect.isclass(obj) and \
                                   not issubclass(obj, Exception))]
all_builtins.append("__import__")                                   

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
    if not hasattr(node, "lineno"): return 0
    return node.lineno
       
#----------------------------------------------------------------------
# Restricted AST nodes & builtins.
#----------------------------------------------------------------------

# Deny evaluation of code if the AST contain any of the following nodes:

unallowed_ast_nodes = [
'AST', 
#'Add', 'And', 'Assert', 'Assign', 'Attribute', 'AugAssign', 'AugLoad', '
#AugStore', 'BinOp', 'BitAnd', 'BitOr', 'BitXor', 
#'BoolOp', 'Break', 'Bytes', 'Call', 'ClassDef', 'Compare', 'Continue', 
#'Del', 'Delete', 'Dict', 'DictComp', 'Div', 'Ellipsis', 'Eq', 
#'ExceptHandler', 'Expr', 'Expression', 'ExtSlice', 'FloorDiv', 'For', 
#'FunctionDef', 'GeneratorExp', 
'Global', 
#'Gt', 'GtE', 'If', 'IfExp', 
'Import', 'ImportFrom', 
#'In', 'Index', 
'Interactive', 
#'Invert', 'Is', 'IsNot', 'LShift', 'Lambda', 'List', 'ListComp', 
#'Load', 
#'Lt', 'LtE', 'Mod', 'Module', 'Mult', 'Name', 
'Nonlocal', 
#'Not', 'NotEq', 'NotIn', 'Num', 'Or', 'Param', 'Pass', 'Pow', 'RShift', 'Raise', 
#'Return', 'Set', 'SetComp', 'Slice', 'Starred', 
#'Store', 'Str', 'Sub', 'Subscript', 'Suite', 'TryExcept', 'TryFinally', 
#'Tuple', 'UAdd', 'USub', 'UnaryOp', 'While', 'With', 'Yield', 
#'alias', 'arg', 'arguments', 'boolop', 'cmpop', 'comprehension', 
#'excepthandler', 'expr', 'expr_context', 'keyword', 'mod', 'operator', 
#'slice', 'stmt', 'unaryop'
]


# Deny evaluation of code if it tries to access any of the following builtins:
unallowed_builtins = [
#'BaseException', 'GeneratorExit', 'KeyboardInterrupt', 'SystemExit', 
#'__build_class__', 
'__import__', 
#'abs', 'all', 'any', 'ascii', 'bin', 'bool', 'bytearray', 'bytes', 
#'chr', 'classmethod', 
'compile', 
#'complex', 
'delattr', 
#'dict', 
'dir', 
#'divmod', 'enumerate', 
'eval', 'exec', 
#'filter', 'float', 'format', 'frozenset', 
'getattr', 'globals', 'hasattr', 
#'hash', 'hex', 'id', 
'input', 
#'int', 'isinstance', 'issubclass', 'iter', 
#'len', 'list', 
'locals', 
#'map', 'max', 'memoryview', 'min', 'next', 'object', 'oct', 
'open', 
#'ord', 'pow', 'print', 'property', 'range', 'repr', 
#'reversed', 'round', 'set', 
'setattr', 
#'slice','sorted', 'staticmethod', 'str', 'sum', 'super', 
#'tuple', 'type', 
'vars', 
#'zip',
]

for ast_name in unallowed_ast_nodes:
    if ast_name == "open": continue #setuptools sandboxing...
    assert(is_valid_ast_node(ast_name))
for name in unallowed_builtins:
    if name == "open": continue #setuptools sandboxing...
    assert(is_valid_builtin(name))

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

class SafeEvalError:
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

class SafeEvalVisitor:
    """
    Data-driven visitor which walks the AST for some code and makes
    sure it doesn't contain any expression/statements which are
    declared as restricted in 'unallowed_ast_nodes'. We'll also make
    sure that there aren't any attempts to access/lookup restricted
    builtin declared in 'unallowed_builtins'. By default we also won't
    allow access to lowlevel stuff which can be used to dynamically
    access non-local envrioments.

    Interface:
      walk(astree) = validate AST and return True if AST is 'safe'

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

    def walk(self, astree):
        "Validate each node in AST and return True if AST is 'safe'."
        self.visit(astree)
        return self.errors == []
        
    def visit(self, node, *args):
        "Recursively validate node and all of its children."
        fn = getattr(self, 'visit' + classname(node))
        if DEBUG: self.trace(node)
        fn(node, *args)
        for child in ast.iter_child_nodes(node):
            self.visit(child, *args)

    def visitName(self, node, *args):
        "Disallow any attempts to access a restricted builtin/attr."
        name = list(ast.iter_fields(node))[0][1]        
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
        for attr in dir(node):
            if attr[:2] != '__':
                print(' ' * 4, "%-15.15s" % attr, getattr(node, attr))

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
      astree = ast.parse(code)
    except SyntaxError:
      astree = None
      pass #will crash anyway later, and we will know the line number
    except:
      print(code)
      raise
    if astree != None:
      checker = SafeEvalVisitor()

      if checker.walk(astree) == False:
        raise SafeEvalCodeException(code, checker.errors)

"""
if __name__ == "__main__":
    s = "a = Integer(3)\n"
    print safe_eval(s)
"""    


