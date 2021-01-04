from sly import Lexer
from errors import ArchimedesError

class ArchimedesLexer(Lexer):
    tokens =  {
        TRUE, FALSE, NIL, MACRO, FLOAT, INT, IDENT, STRING,
        EXP, DIV, MUL, ADD, SUB, ASSIGN, GT, GTE, LT,
        LTE, MOD, VAR, OR, AND, EQ, NEQ, DOLLAR,
        QUESTION, MATCH, EXC, RPAREN, LPAREN, RBRACKET,
        LBRACKET, COMMA, COLON, BIND, FUNC
    }

    def error(self, t):
        raise ArchimedesError('Syntax error in line %d, near token %s' % (self.lineno, t.value))

    @_(r'\n+')
    def ignore_newline(self, t):
        self.lineno += t.value.count('\n')

    ignore_space = r"(\s){1,}"
    ignore_comment = r'#.*$'

    @_(r'[0-9]+\.[0-9]+(e[+-]?[0-9]+)?')
    def FLOAT(self, t):
         t.value = float(t.value)
         return t

    @_(r'[0-9]+')
    def INT(self, t):
         t.value = int(t.value)
         return t

    @_(r"'(?:[^']|'')*'",
       r'"(?:[^"]|"")*"')
    def STRING(self, t):
        if t.value[0] == "'":
            t.value = t.value[1:-1].replace(r"''", "'")
        else:
            t.value = t.value[1:-1].replace(r'""', '"')
        return t

    TRUE = r'true'
    FALSE = r'false'
    NIL = r'nil'
    IDENT = r'[A-Za-z][\.\w]*'
    EXP = r'\^'
    DIV = r'\/'
    MUL = r'\*'
    ADD = r'\+'
    SUB = r'\-'
    FUNC = r'=>'
    GTE = r'>='
    GT = r'>'
    LTE = r'<='
    LT = r'<'
    MOD = r'\%'
    VAR = r'\@'
    OR = r'\|\|'
    AND = r'\&\&'
    EQ = r'=='
    NEQ = r'!='
    DOLLAR = r'\$'
    QUESTION = r'\?'
    MATCH = r'=\~'
    ASSIGN = r'='
    EXC = r'\!'

    RPAREN = r'\)'
    LPAREN = r'\('
    RBRACKET = r'\]'
    LBRACKET = r'\['
    COMMA = r'\,'
    BIND = r'\:\:'
    COLON = r'\:'


    #rule(/\{(.*?)\}/m) { |t| [ :MACRO, t[1..-2] ] }
