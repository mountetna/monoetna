from typing import Any
from collections.abc import Iterable

def _ensure_base(val: Any, fxn):
    # Treat strings and bytes as non-iterables in this context
    if isinstance(val, Iterable) and not isinstance(val, (str, bytes)):
        return fxn(val)
    else:
        return fxn([val])

def _all(val: Any):
    """
    Implementation of the base 'all()' function that is safe to use with non-iterables
    """
    return _ensure_base(val, all)

def _any(val: Any):
    """
    Implementation of the base 'any()' function that is safe to use with non-iterables
    """
    return _ensure_base(val, any)