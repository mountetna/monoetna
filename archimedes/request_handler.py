from errors import ArchimedesError
from parse import ArchimedesParser
from lex import ArchimedesLexer

def execute(request):
    """
    Executing requests using the methods defined in the methods module.
    
    Request dict format:
    {"func" : <function name that is in the methods module>,
     "requests" : <list of dicts, functions and arguments to run in order. Not used if "func" is not equal to "tether">,
     "args" : <ordered arguments for func, including data input. For tether, these are the initial input arguments>,
     "kwargs" : <key worded arguments for func>
    }
    
    input: request (dict)
    output: <tether function> OR method.func(*args, **kwargs)
    """

    token = request.headers.get('Etna-Token')
    params = request.get_json(force=True)

    project_name = params['project_name']
    manifest = params['manifest']
    result = resolve(manifest, token)
    return { "result": result }

def resolve(query, token=None):
    arch_lexer = ArchimedesLexer()
    arch_parser = ArchimedesParser()
    #import code; code.interact(local=dict(globals(), **locals()))
    arch_parser.parse(arch_lexer.tokenize(query))
    return arch_parser.return_vars
