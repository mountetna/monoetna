from typing import Callable, Mapping, Any

import inspect


def inject(fn: Callable, injections: Mapping[str, Any]):
    sig: inspect.Signature = inspect.signature(fn)
    parameters = sig.parameters.values()
    inject_p = {}

    for param in parameters:
        if param.name in injections:
            inject_p[param.name] = injections[param.name]
        if param.kind == param.VAR_KEYWORD:
            inject_p.update(injections)
            break

    return fn(**inject_p)
