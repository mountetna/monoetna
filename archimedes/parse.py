from sly import Parser
from lex import ArchimedesLexer
from pandas import Series, DataFrame, isna, concat
from errors import ArchimedesError
from vector import Vector
from matrix import Matrix
import re

class Op:
    def __init__(self, *args, **opts):
        self.args = args
        self.func = opts.get('func')
        self.name = opts.get('name')

    def resolve(self):
        return self.func(*self.args) if self.func else self.args[0]

class ArchimedesParser(Parser):
    tokens = ArchimedesLexer.tokens
    debugfile = 'parser.out'

    def error(self,p):
        if p:
             msg = "Syntax error in line %d, near %s" % (p.lineno, p.value)
        else:
             msg = "Syntax error at EOF"

        raise ArchimedesError(msg)

    def __init__(self):
        self.names = {}
        self.vars = {}
        self.return_vars = {}

    precedence = (
        ('left', QUESTION),
        ('left', MOD),
        ('left', GT, GTE, LT, LTE, EQ, MATCH),
        ('left', SUB, ADD),
        ('left', DIV, MUL),
        ('left', EXP),
        ('left', DOLLAR),
        ('left', VAR)
    )

    start = 'script'

    @_('')
    def empty(self,p):
        pass

    @_('assignment')
    def script(self, p):
        return Op(p.assignment, func=lambda a : [a()], name='script').resolve

    @_('assignment script')
    def script(self, p):
        return Op(p.assignment, p.script, func=lambda a, s : [a()] + s(), name='script').resolve

    @_('VAR IDENT ASSIGN e')
    def assignment(self, p):
        def assign(i,e,v,r):
            e = e()
            v[i] = e
            r[i] = e
            return e
        return Op(p.IDENT, p.e, self.vars, self.return_vars, func=assign, name='assignment').resolve

    @_('IDENT ASSIGN e')
    def assignment(self, p):
        def assign(i,e,v):
            v[i] = e()
            return v[i]
        return Op(p.IDENT, p.e, self.vars, func=assign, name='assignment').resolve

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

    @_('IDENT', 'VAR IDENT')
    def var(self, p):
        def getvar(i,v,l):
            if not i in self.vars:
                raise ArchimedesError("Syntax error in line %d: No such variable %s" % ( l, i) )
            return v[i]

        return Op(p.IDENT, self.vars, p.lineno, func=getvar).resolve

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

    @_('LPAREN func_args RPAREN FUNC LBRACKET func_block RBRACKET')
    def e(self,p):
        return Op(p.func_args, p.func_block, func=lambda a, b : Function(a, b)).resolve

    @_('func_arg COMMA func_args')
    def func_args(self, p):
        return [ p.func_arg ] + p.func_args

    @_('empty')
    def func_args(self, p):
        return []

    @_('empty')
    def func_block(self, p):
        return []

    @_('IDENT')
    def func_arg(self, p):
        return FuncArg( p.IDENT )

    @_('IDENT EQ e')
    def func_arg(self, p):
        return FuncArg(p.IDENT, e)

    @_('var LPAREN args RPAREN')
    def e(self, p):
        return self.function(p.var, p.args)

    @_('var LPAREN RPAREN')
    def e(self, p):
        return self.function(p.var, [])

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
        return Op(p.e, func=lambda e : e())

    @_('empty')
    def slice(self,p):
        return Op(None)
