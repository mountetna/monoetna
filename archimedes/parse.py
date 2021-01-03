from sly import Parser
from lex import ArchimedesLexer
from pandas import Series, DataFrame, isna, concat
from errors import ArchimedesError
from vector import Vector
from matrix import Matrix
import re

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
        return [p.assignment]

    @_('assignment script')
    def script(self, p):
        return [ p.assignment ] + p.script

    @_('VAR IDENT ASSIGN e')
    def assignment(self, p):
        self.vars[p.IDENT] = p.e
        self.return_vars[p.IDENT] = p.e

    @_('IDENT ASSIGN e')
    def assignment(self, p):
        self.vars[p.IDENT] = p.e

    @_('FLOAT')
    def e(self, p):
        return p.FLOAT

    @_('INT')
    def e(self, p):
        return p.INT

    @_('STRING')
    def e(self, p):
        return p.STRING

    @_('vector')
    def e(self, p):
        return p.vector

    # some defined keywords
    @_('TRUE')
    def e(self, p):
        return True
    @_('FALSE')
    def e(self, p):
        return False

    @_('NIL')
    def e(self, p):
        return None

    @_('IDENT', 'VAR IDENT')
    def var(self, p):
        if not p.IDENT in self.vars:
            raise ArchimedesError("Syntax error in line %d: No such variable %s" % ( p.lineno, p.IDENT ) )

        return self.vars[p.IDENT]

    @_('var')
    def e(self, p):
        return p.var

    @_('EXC e')
    def e(self, p):
        return not p.e

    # Basic math operations, including the ternary operator
    @_('LPAREN e RPAREN')
    def e(self, p):
        return p.e

    @_('e ADD e')
    def e(self, p):
        return p.e0 + p.e1

    @_('e EXP e')
    def e(self, p):
        return p.e0 ** p.e1

    @_('e SUB e')
    def e(self, p):
        return p.e0 - p.e1

    @_('e QUESTION e COLON e')
    def e(self, p):
        if isinstance(p.e0, Series):
            return p.e0.ternary(p.e1,p.e2)
        else:
            return p.e1 if p.e0 else p.e2

    # indexing, used by Vectors and Matrices
    @_('e LBRACKET e RBRACKET')
    def e(self, p):
        return p.e0[p.e1]

    # matrix slicing
    @_('e LBRACKET slice COMMA slice RBRACKET')
    def e(self, p):
        return p.matrix.slice(rows, cols)

    # Comparison operators
    @_('e GT e')
    def e(self, p):
        return p.e0 > p.e1

    @_('e GTE e')
    def e(self, p):
        return p.e0 >= p.e1

    @_('e LT e')
    def e(self, p):
        return p.e0 < p.e1

    @_('e LTE e')
    def e(self, p):
        return p.e0 <= p.e1


    # Matrix column reference notation
    @_('e DOLLAR IDENT')
    def e(self, p):
        return p.e[p.IDENT]


    # Arithmetic operations
    @_('SUB e')
    def e(self, p):
        return -p.e

    @_('e MOD e')
    def e(self, p):
        return p.e0 % p.e1

    @_('e DIV e')
    def e(self, p):
        return p.e0 / p.e1

    @_('e MUL e')
    def e(self, p):
        return p.e0 * p.e1

    @_('e OR e')
    def e(self, p):
        if (isinstance(p.e0, Series) or isinstance(p.e1, Series)):
            return p.e0 | p.e1
        return p.e0 or p.e1

    @_('e AND e')
    def e(self, p):
        if (isinstance(p.e0, Series) or isinstance(p.e1, Series)):
            return p.e0 & p.e1
        return p.e0 and p.e1

    @_('e EQ e')
    def e(self, p):
        if isinstance(p.e0, Series) and p.e1 is None:
            return isna(p.e0)
        if isinstance(p.e1, Series) and p.e0 is None:
            return isna(p.e1)
        return p.e0 == p.e1

    @_('e NEQ e')
    def e(self, p):
        if isinstance(p.e0, Series) and p.e1 is None:
            return ~isna(p.e0)
        if isinstance(p.e1, Series) and p.e0 is None:
            return ~isna(p.e1)
        return p.e0 != p.e1

    @_('e MATCH e')
    def e(self, p):
        return re.match(p.e1, p.e0)

    @_('BIND e')
    def e(self, p):
        return Matrix(p.e)

    @_('e BIND e')
    def e(self, p):
        if isinstance(p.e0, DataFrame):
            e1 = Matrix(p.e1)
            e1.index = p.e0.index
            return concat([ p.e0, e1 ], axis = 1)
        raise ArchimedesError('You must bind to an existing matrix')

    @_('var LPAREN args RPAREN')
    def e(self, p):
        return self.function(p.var, p.args)

    @_('var LPAREN RPAREN')
    def e(self, p):
        return self.function(p.var, [])

    @_('LBRACKET vector_items RBRACKET')
    def vector(self, p):
        return Series(
            data=[ data for (index,data) in p.vector_items ],
            index=[ index for (index,data) in p.vector_items ]
        )

    @_('LBRACKET empty RBRACKET')
    def vector(self, p):
        return Series([])

    @_('vector_item')
    def vector_items(self, p):
        return [ p.vector_item ]

    @_('vector_item COMMA vector_items')
    def vector_items(self, p):
        return [ p.vector_item ] + p.vector_items

    @_('e')
    def vector_item(self, p):
        return ( None, p.e )

    @_('IDENT COLON e')
    def vector_item(self, p):
        return ( p.IDENT, p.e )

    @_('STRING COLON e')
    def vector_item(self, p):
        return ( p.STRING, p.e )

    @_('e')
    def args(self, p):
        return [p.e]

    @_('e COMMA args')
    def args(self, p):
        return [p.e] + p.args

    @_('e')
    def slice(self,p):
        return p.e

    @_('empty')
    def slice(self,p):
        return None
