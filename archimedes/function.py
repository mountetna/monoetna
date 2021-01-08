from errors import ArchimedesError
from scope import Scope

class Function:
    def __init__(self, parser, arg_defs, production):
        self.arg_defs = arg_defs
        self.production = production
        self.parent_scope = parser.current_scope
        self.parser = parser

    def getargs(self, args):
        computed_args = []
        for i, (name, default) in enumerate(self.arg_defs):
            if i < len(args):
                computed_args.append( (name, args[i]) )
            elif default:
                computed_args.append( (name, default()) )
            else:
                raise ArchimedesError("Missing argument %s with no default" % (name))

        return computed_args
            

    def call(self, *args):
        # create a new scope
        previous_scope = self.parser.current_scope

        scope = Scope(self.parent_scope)

        self.parser.current_scope = scope

        # set the arguments
        for (arg_var, value) in self.getargs(args):
            scope.set(arg_var, value)

        # call the script
        value = self.production()

        # set the current scope back to the previous scope
        self.parser.current_scope = previous_scope

        return value
