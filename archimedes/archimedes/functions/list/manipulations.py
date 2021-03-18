import itertools

def unique(ls):
    return list(dict.fromkeys(ls))

def flatten(ls):
    return list(itertools.chain.from_iterable(ls))