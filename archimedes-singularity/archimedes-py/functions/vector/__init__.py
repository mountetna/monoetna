from pandas import Series
from numpy import arange, unique


def length(vec):
    return vec.size


def max(vec):
    return vec.max()


def min(vec):
    return vec.min()


def rep(vector, times):
    return vector.repeat(times)


def join(vector, sep=""):
    return sep.join(vector.astype(str))


def seq(start, stop, interval=1):
    return Series(arange(start, stop + 1, interval))


def unique(vec):
    return Series(vec.unique())


def compact(vec):
    return vec.dropna()


def filter(vec, func):
    return Series(
        [vec[i] for i in range(vec.size) if func.call(vec[i], i, vec.index[i])]
    )


def any(vec, func):
    return bool(vec.apply(lambda x: func.call(x)).any())


def all(vec, func):
    return bool(vec.apply(lambda x: func.call(x)).all())


def group(vector, func):
    groups = {}
    for v in vector:
        group = func.call(v)
        if not isinstance(group, str):
            raise ValueError("Must group by string factors!")
        if not group in groups:
            groups[group] = []
        groups[group].append(v)

    return Series([Series(v) for v in groups.values()], index=groups.keys())


def set_labels(vec, labels):
    vec = vec.copy()
    vec.index = labels
    return vec


def labels(vec):
    return Series(vec.index.tolist())
