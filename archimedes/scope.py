class Scope:
    def __init__(self, parent_scope=None):
        self.parent_scope = parent_scope
        self.vars = {}
        self.export_vars = {}

    def set(self, var_name, value):
        self.vars[var_name] = value
        return value

    def export_set(self, var_name, value):
        self.export_vars[var_name] = value

        return value

    def export_get(self, var_name):
        return self.export_vars[var_name]

    def get(self, var_name):
        if not self.parent_scope or var_name in self.vars:
            return self.vars[var_name]

        return self.parent_scope.get(var_name)
