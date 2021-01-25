from pylint.checkers import BaseChecker
from pylint.interfaces import IAstroidChecker
from pylint.lint import PyLinter
from pylint.reporters.text import TextReporter


class ArchimedesDslChecker(BaseChecker):
    __implements__ = IAstroidChecker

    name = "archmedes-dsl"
    priority = -100
    msgs = {
        "W9001": (
            "Use only safe imports.",
            "non-safe-imports",
            "Scripts should only import from archimedes.functions",
        ),
        "W9002": (
            "Use no relative imports",
            "no-relative-imports",
            "Scripts should not use any relative imports",
        ),
        "W9003": (
            "Use no unsafe globals or unsafe names",
            "no-unsafe-globals",
            "Scripts should not use any unsafe globals",
        ),
        "W9004": (
            "Use no private members",
            "no-private-members",
            "Scripts should not use any private members",
        ),
    }
    options = tuple()

    def visit_import(self, node):
        for (module, as_name) in node.names:
            if not module.startswith("archimedes.functions."):
                self.add_message("non-safe-imports", node=node)
                break

    def visit_importfrom(self, node):
        if node.level:
            self.add_message("no-relative-imports", node=node)
            return

        if not node.modname.startswith("archimedes.functions."):
            self.add_message("non-safe-imports", node=node)

    # NOTE: This is NOT the security model.  A clever attack can still access the full
    # runtime with clever usage of builtins.  However, restricting access to these identifiers
    # encourages users to stick to our intended library API so that we don't end up
    # having assumptions being built on our environment or using tools we don't want to
    # support.  The security model is handled by the runner module, which integrates with
    # docker processes to expose a sandbox that is limited in its capabilities.
    UNSAFE_NAMES = {
        "getattr",
        "setattr",
        "delattr",
        "hasattr",
        "compile",
        "__builtins__",
        "__module__",
        "__name__",
        "eval",
        "exec",
        "locals",
        "globals",
        "__package__",
        "__spec__",
        "__loader__",
        "__import__",
        "__debug__",
        "__doc__",
        "__build_class__",
        "exit",
        "quit",
    }

    def visit_name(self, node):
        if node.name in self.UNSAFE_NAMES:
            self.add_message("no-unsafe-globals", node=node)

    def visit_attribute(self, node):
        if node.attrname.startswith("_"):
            self.add_message("no-private-members", node=node)


def run(files, reporter=None):
    if not reporter:
        reporter = TextReporter()
        reporter.set_output()

    old_writeln = reporter.writeln
    # Shared state.  Python closures are -final-, that is not mutable by reference.
    # However, mutating inner state of the final reference is allowed, so we use a
    # 1 length array to capture this.
    success = [True]
    # If you can believe it, pylint doesn't actually give you any data back on wether a linting
    # run is successful, it expects you to parse output to determine that :'(
    # To compensate, we monkey patch writeln, and flag a failure as any linter resulting in
    # output.
    def new_writeln(*args, **kwds):
        success[0] = False
        return old_writeln(*args, **kwds)

    reporter.__setattr__("writeln", new_writeln)

    linter = PyLinter(reporter=reporter)
    linter.register_checker(ArchimedesDslChecker(linter))
    linter.check(files)
    return success[0]


def main():
    import sys

    if not run([sys.argv[1]]):
        sys.exit(1)


if __name__ == "__main__":
    main()
