from collections.abc import Iterable
from typing import List, Union
import re 

def ensure_list(l):
    if not isinstance(l, list):
        return [l]
    return l

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

def reverse(ls: list):
    return [element for element in reversed(ls)]

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

def order(ls: Union[List[Union[int, float]], List[str], List[bool]], return_indexes: bool = False, decreasing: bool = False):
    """
    Output: A list. ls reordered in "increasing" order unless 'decreasing' is set to "False"
    types: str, bool, 'number'
    """
    if all(map(lambda x: isinstance(x, str), ls)):
        output=order_str(ls, return_indexes)
    elif all(map(lambda x: isinstance(x, bool), ls)):
        num_ls = list({True: 0, False: 1}[x] for x in ls)
        order_inds = order_num(num_ls, True)
        if return_indexes: output=order_inds
        else: output=list(ls[order_inds])
    elif all(map(lambda x: isinstance(x, (int, float)), ls)):
        output=order_num(ls, return_indexes)
    elif all(map(lambda x: isinstance(x, type(ls[0])), ls)):
        raise NotImplemented("Our 'order' function is not yet built for your data type. *Please let the Data Library team know how you got here!*")
    else:
        raise Exception("Cannot 'order' mixed data types. Can the data be cleaned?")
    return reverse(output) if decreasing else output
