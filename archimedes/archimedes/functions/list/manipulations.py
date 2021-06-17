from collections.abc import Iterable
from typing import List, Union
import re 

def unique(ls: list):
    return list(dict.fromkeys(ls))

def flatten(ls: list):
    def _flatten(l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from _flatten(el)
            else:
                yield el
    return list(_flatten(ls))

def order_num(ls: List[Union[int, float]], return_indexes: bool = False):
    if return_indexes:
        ret = 0
    else:
        ret = 1
    return list(i[ret] for i in sorted(enumerate(ls), key = lambda x:x[1]))

def order_str(ls: List[str], return_indexes: bool = False): 
    """ Sort a string iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key[1].lower()) ]
    if return_indexes:
        ret = 0
    else:
        ret = 1
    return list(i[ret] for i in sorted(enumerate(ls), key = alphanum_key))

def order(ls: list, return_indexes: bool = False):
    """
    Output: ls reordered in "increasing" order
    types: str, bool, 'number' (else case):
    """
    if all(map(lambda x: isinstance(x, str), ls)):
        return order_str(ls, return_indexes)
    elif all(map(lambda x: isinstance(x, bool), ls)):
        num_ls = list({True: 0, False: 1}[x] for x in ls)
        order_inds = order_num(num_ls, True)
        if return_indexes: return order_inds
        else: return ls[order_inds]
    else:
        return order_num(ls, return_indexes)
