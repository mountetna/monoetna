from pandas import DataFrame

def Matrix(s):
    return DataFrame(
        { label or idx : column.tolist() for (idx, (label,column)) in enumerate(s.iteritems()) },
        index = s[0].index
    )
