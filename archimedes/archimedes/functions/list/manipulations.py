from collections.abc import Iterable

def unique(ls):
    return list(dict.fromkeys(ls))

# def flatten(ls):
#     return list(itertools.chain.from_iterable(ls))

def flatten(ls):
    def _flatten(l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el
    return list(_flatten(ls))