from collections.abc import Iterable
import re 

def unique(ls):
    return list(dict.fromkeys(ls))

def flatten(ls):
    def _flatten(l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el
    return list(_flatten(ls))

def order(ls): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key.lower()) ] 
    return sorted(ls, key = alphanum_key)