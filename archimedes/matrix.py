from pandas import DataFrame

def Matrix(s):
    return DataFrame(
        { label : column.tolist() for (label,column) in s.iteritems() },
        index = s[0].index
    )
