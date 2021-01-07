from errors import ArchimedesError

class Function:
    def __init__(self, arg_defs, production, parent_scope):
        self.arg_defs = arg_defs
        self.production = production
        self.parent_scope = parent_scope

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
            

    def call(self):
        return self.production()
