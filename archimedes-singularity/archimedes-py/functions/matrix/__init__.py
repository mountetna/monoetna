from pandas import Series
import pandas as pd


def spread(matrix):
    return matrix


def rowsums(matrix):
    return matrix.sum(0)


def colsums(matrix):
    return matrix.sum(1)


def rownames(matrix):
    return Series(matrix.index.tolist())


def colnames(matrix):
    return Series(matrix.columns.tolist())


def numrows(matrix):
    return matrix.index.size


def numcols(matrix):
    return matrix.columns.size
