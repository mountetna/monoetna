from sly import Parser
from lex import ArchimedesLexer
from pandas import Series, DataFrame, isna, concat
from errors import ArchimedesError
from vector import Vector
from matrix import Matrix
from function import Function
from scope import Scope
import inspect
import re
import sys

class Op:
    def __init__(self, *args, **opts):
        self.args = args
        self.func = opts.get('func')

        #caller = inspect.stack()[1]
        #self.name = caller[3]
        #self.code = caller[4][0]

    def resolve(self):
        return self.func(*self.args) if self.func else self.args[0]


class ArchimedesParser(Parser):
    tokens = ArchimedesLexer.tokens
    debugfile = 'parser.out'

    def payload(self):
        return self.root_scope.export_vars

    def error(self,p):
        if p:
             msg = "Syntax error in line %d, near %s" % (p.lineno, p.value)
        else:
             msg = "Syntax error at EOF"

        raise ArchimedesError(msg)

    def __init__(self):
        self.names = {}
        self.root_scope = Scope()

        self.current_scope = self.root_scope
        self.return_var = None

    precedence = (
        ('right',FUNC),
        ('right', ASSIGN),
        ('left', QUESTION),
        ('left', MOD),
        ('left', GT, GTE, LT, LTE, EQ, MATCH),
        ('left', SUB, ADD),
        ('left', DIV, MUL),
        ('left', EXP),
        ('left', DOLLAR),
        ('left', VAR),
    )

    start = 'script'

    @_('')
    def empty(self,p):
        pass

    @_('script_items')
    def script(self,p):
        def run_script(script_items):
            for script_item in script_items:
                # Assignment or return
                script_item()

                if self.return_var:
                    if self.current_scope == self.root_scope:
                        raise ArchimedesError("Cannot return from root scope!")

                    value = self.return_var[0]
                    self.return_var = None
                    return value

        return Op(p.script_items, func=run_script).resolve

    @_('script_item')
    def script_items(self, p):
        return [ p.script_item ]

    @_('script_item script_items')
    def script_items(self, p):
        return [ p.script_item ] + p.script_items

    @_('assignment')
    def script_item(self, p):
        return Op(p.assignment, func=lambda a : a()).resolve

    @_('RETURN e')
    def script_item(self,p):
        def returns(e):
            self.return_var = [ e() ]
        return Op(p.e, func=returns).resolve

    @_('VAR IDENT ASSIGN e')
    def assignment(self, p):
        def export_assign(ident,e,lineno):
            if self.current_scope != self.root_scope:
                raise ArchimedesError("Syntax error in line %d: Cannot export variables from non-root scope" % (l))
            return self.current_scope.export_set(ident,e())
        return Op(p.IDENT, p.e, p.lineno, func=export_assign).resolve

    @_('IDENT ASSIGN e')
    def assignment(self, p):
        def assign(i,e):
            return self.current_scope.set(i,e())

        return Op(p.IDENT, p.e, func=assign).resolve

    @_('FLOAT')
    def e(self, p):
        return Op(p.FLOAT).resolve

    @_('INT')
    def e(self, p):
        return Op(p.INT).resolve

    @_('STRING')
    def e(self, p):
        return Op(p.STRING).resolve

    @_('vector')
    def e(self, p):
        return Op(p.vector, func=lambda v: v()).resolve

    # some defined keywords
    @_('TRUE')
    def e(self, p):
        return Op(True).resolve

    @_('FALSE')
    def e(self, p):
        return Op(False).resolve

    @_('NIL')
    def e(self, p):
        return Op(None).resolve

    @_('LPAREN arg_defs RPAREN FUNC LBRACE script RBRACE')
    def e(self,p):
        return Op(p.arg_defs, p.script, func=lambda arg_defs, script : Function(arg_defs, script, self.current_scope)).resolve

    @_('LPAREN arg_defs RPAREN FUNC e')
    def e(self,p):
        return Op(p.arg_defs, p.e, func=lambda arg_defs, e : Function(arg_defs, e, self.current_scope)).resolve

    @_('arg_def')
    def arg_defs(self, p):
        return [ p.arg_def ]

    @_('arg_def COMMA arg_defs')
    def arg_defs(self, p):
        return [ p.arg_def ] + p.arg_defs

    @_('empty')
    def arg_defs(self, p):
        return []

    @_('IDENT')
    def arg_def(self, p):
        return ( p.IDENT, None )

    @_('IDENT ASSIGN e')
    def arg_def(self, p):
        return ( p.IDENT, p.e )

    @_('IDENT')
    def var(self, p):
        def getvar(i,l):
            try:
                return self.current_scope.get(i)
            except KeyError as e:
                raise ArchimedesError("Syntax error in line %d: No such variable %s" % ( l, i) )
        return Op(p.IDENT, p.lineno, func=getvar).resolve

    @_('VAR IDENT')
    def var(self, p):
        def getvar(i,l):
            try:
                return self.root_scope.export_get(i)
            except KeyError as e:
                raise ArchimedesError("Syntax error in line %d: No such exported variable %s" % ( l, i) )

        return Op(p.IDENT, p.lineno, func=getvar).resolve

    @_('var')
    def e(self, p):
        return Op(p.var, func=lambda v : v()).resolve

    @_('EXC e')
    def e(self, p):
        return Op(p.e, func=lambda e : not e()).resolve

    # Basic math operations, including the ternary operator
    @_('LPAREN e RPAREN')
    def e(self, p):
        return Op(p.e, func=lambda e : e()).resolve

    @_('e ADD e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() + e1()).resolve

    @_('e EXP e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() ** e1()).resolve

    @_('e SUB e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() - e1()).resolve

    @_('e QUESTION e COLON e')
    def e(self, p):

        def ternary(e0, e1, e2):
            e0 = e0()
            if isinstance(e0, Series):
                return None
            else:
                return e1() if e0 else e2()

        return Op(p.e0, p.e1, p.e2,func=ternary).resolve

    # indexing, used by Vectors and Matrices
    @_('e LBRACKET e RBRACKET')
    def e(self, p):
        return Op(p.e0,p.e1, func=lambda e0, e1 : e0()[e1()]).resolve

    # matrix slicing
    @_('e LBRACKET slice COMMA slice RBRACKET')
    def e(self, p):
        return Op(p.e, p.slice1, p.slice2, func=lambda e, s1, s2 : e()[s1(), s2()]).resolve

    # Comparison operators
    @_('e GT e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() > e1()).resolve

    @_('e GTE e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() >= e1()).resolve

    @_('e LT e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() < e1()).resolve

    @_('e LTE e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() <= e1()).resolve


    # Matrix column reference notation
    @_('e DOLLAR IDENT')
    def e(self, p):
        return Op(p.e, p.IDENT, func=lambda e, i : e()[i]).resolve

    # Arithmetic operations
    @_('SUB e')
    def e(self, p):
        return Op(p.e, func=lambda e : -e()).resolve

    @_('e MOD e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() % e1()).resolve

    @_('e DIV e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() / e1()).resolve

    @_('e MUL e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : e0() * e1()).resolve

    @_('e OR e')
    def e(self, p):
        def orf(e0, e1):
            e0 = e0()
            e1 = e1()
            if (isinstance(e0, Series) or isinstance(e1, Series)):
                return e0 | e1
            else:
                return e0 or e1
        return Op(p.e0, p.e1, func=orf).resolve

    @_('e AND e')
    def e(self, p):
        def andf(e0, e1):
            e0 = e0()
            e1 = e1()
            if (isinstance(e0, Series) or isinstance(e1, Series)):
                return e0 & e1
            else:
                return e0 and e1
        return Op(p.e0, p.e1, func=andf).resolve

    @_('e EQ e')
    def e(self, p):
        def equals(e0, e1):
            e0 = e0()
            e1 = e1()
            if isinstance(e0, Series) and e1 is None:
                return isna(e0)
            if isinstance(e1, Series) and e0 is None:
                return isna(e1)
            return e0 == e1
        return Op(p.e0, p.e1, func=equals).resolve

    @_('e NEQ e')
    def e(self, p):
        def nequals(e0, e1):
            e0 = e0()
            e1 = e1()
            if isinstance(e0, Series) and e1 is None:
                return ~isna(e0)
            if isinstance(e1, Series) and e0 is None:
                return ~isna(e1)
            return e0 != e1
        return Op(p.e0, p.e1, func=nequals).resolve

    @_('e MATCH e')
    def e(self, p):
        return Op(p.e0, p.e1, func=lambda e0, e1 : re.match(e1(), e0())).resolve

    @_('BIND e')
    def e(self, p):
        return Op(p.e, func=lambda e : Matrix(e())).resolve

    @_('e BIND e')
    def e(self, p):
        def bind(e0, e1):
            e0 = e0()
            e1 = e1()
            if isinstance(e0, DataFrame):
                e1 = Matrix(e1)
                e1.index = e0.index
                return concat([ e0, e1 ], axis = 1)
            raise ArchimedesError('You must bind to an existing matrix')
        return Op(p.e0, p.e1, func=bind).resolve

    @_('e LPAREN args RPAREN')
    def function_call(self, p):
        return ( p.e, p.args )

    @_('e LPAREN RPAREN')
    def function_call(self, p):
        return ( p.e, None )

    @_('function_call')
    def e(self, p):
        def call_function(func, args):
            func = func()
            if args:
                args = args()
            else:
                args = []
            if not isinstance(func, Function):
                raise ArchimedesError("Syntax error: not a function" % (l))

            # create a new scope
            previous_scope = self.current_scope

            self.current_scope = Scope(func.parent_scope)

            # set the arguments
            for (arg_var, value) in func.getargs(args):
                self.current_scope.set(arg_var, value)

            # call the script
            value = func.call()

            # set the current scope back to the previous scope
            self.current_scope = previous_scope

            return value

        return Op(*p.function_call, func=call_function).resolve

    @_('LBRACKET vector_items RBRACKET')
    def vector(self, p):
        def vec(vi):
            vi = vi()
            return Series(
                data=[ data for (index,data) in vi ],
                index=[ index for (index,data) in vi ]
            )
        return Op(p.vector_items, func=vec).resolve

    @_('LBRACKET empty RBRACKET')
    def vector(self, p):
        return Op(None, func=lambda : Series([])).resolve

    @_('vector_item')
    def vector_items(self, p):
        return Op(p.vector_item, func=lambda vi : [ vi() ]).resolve

    @_('vector_item COMMA vector_items')
    def vector_items(self, p):
        return Op(p.vector_item, p.vector_items, func=lambda vi, vis : [ vi() ] + vis()).resolve

    @_('e')
    def vector_item(self, p):
        return Op(p.e, func=lambda e : ( None, e() )).resolve

    @_('IDENT COLON e')
    def vector_item(self, p):
        return Op( p.IDENT, p.e, func=lambda i, e : ( i, e() ) ).resolve

    @_('STRING COLON e')
    def vector_item(self, p):
        return Op( p.STRING, p.e, func=lambda s, e : ( s, e() ) ).resolve

    @_('e')
    def args(self, p):
        return Op( p.e, func=lambda e : [e()] ).resolve

    @_('e COMMA args')
    def args(self, p):
        return Op(p.e, p.args, func=lambda e, a : [e()] + args()).resolve

    @_('e')
    def slice(self,p):
        return Op(p.e, func=lambda e : e()).resolve

    @_('empty')
    def slice(self,p):
        return Op(None).resolve
