import numpy as np

DISCRETE_KINDS = 'ObUS'
CONTINUOUS_KINDS = 'ifuc'

def is_discrete(nda: np.ndarray):
    return nda.dtype.kind in DISCRETE_KINDS

def is_continuous(nda: np.ndarray):
    return nda.dtype.kind in CONTINUOUS_KINDS

def scale(y):
    x = y.copy()
    x -= x.mean()
    x /= x.std()
    return x

def _leave_default_or_null(target, default, none_if = False, default_when = "make"):
    if target == default_when:
        if none_if:
            return None
        else:
            return default
    return target

def _isCol(test, data_frame, return_values=False):
    '''
    This function tests if an input is the name of a meta.data slot in a target data_frame.
    A carryover from dittoSeq where the check was of a sub-slot of a different part of multiple object types.
    Might remove, but keeping for now while I port from R to python.
    
    'test' String or vector of strings, the "potential.metadata.name"(s) to check for.
    'data_frame' A data.frame.
    'return_values' Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}

    Value: Returns list of bool(s) indicating whether each instance in 'test' is a column of the 'data_frame'.
    Alternatively, returns the values of 'test' that were indeed columns of 'data_frame' if 'return_values = True'.

    Examples
    .isCol("timepoint", df) # True

    # To test if many things are metadata of an data_frame
    .isCol(["age","groups"], df) # [False, True]

    # 'return_values' input is especially useful in these cases.
    .isCol(["age","groups"], df, return_values = True) # "groups"
    '''
    if (return_values):
        return [t for t in test if _isCol(t, data_frame, return_values=False)[0]]
    cols = data_frame.columns
    try:
        iterator = iter(test)
    except TypeError:
        return list(test in cols)
    else:
        return list(x in cols for x in test)

def _col(col, data_frame, adjustment = None, adj_fxn = None):
    '''
    Returns the values of a column, where values are named by the rownames of data_frame

    'col' String, the name of the column to grab.
    'data_frame' A data.frame.
    'adjustment' Ignored if the target metadata is not numeric, a recognized string indicating whether numeric data should be used directly (default) versus adjusted to be
        "z-score": scaled with the scale() function to produce a relative-to-mean z-score representation
        "relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]
    'adj.fxn' A lambda function which takes an array of values and returns an array of the same length.
        For example, `lambda x: log2(x)`.

    Value: A list.
    @details
    Retrieves the values of a metadata slot from \code{object}, or the clustering slot if \code{meta = "ident"} and the \code{object} is a Seurat.

    If \code{adjustment} or \code{adj.fxn} are provided, then these requested adjustments are applied to these values (\code{adjustment} first).
    Note: Alterations via \code{adjustment} are only applied when metadata is numeric, but \code{adj.fxn} alterations are applied to metadata of any type.

    Lastly, outputs these values are named as the cells'/samples' names.
    
    examples:
    example(importDittoBulk, echo = FALSE)
    ._col("groups", object = myRNA)

    myRNA$numbers <- seq_len(ncol(myRNA))
    ._col("numbers", myRNA, adjustment = "z-score")
    ._col("numbers", myRNA, adj.fxn = as.factor)
    ._col("numbers", myRNA, adj.fxn = function(x) \{log2(x)\})
    '''
    if not _isCol(col, data_frame):
        raise Exception("" + col + " is not a column of 'data_frame'.")
    # Retrieve target columns's values and make into numpy array
    values = np.array(data_frame[col])
    # Add 'adjustment's, if dtype is numeric.
    if (adjustment != None and is_continuous(values)):
        if (adjustment=="z-score"):
             values = scale(values)
        if (adjustment=="relative.to.max"):
             values = values/max(values)
    if adj_fxn!=None:
        values = adj_fxn(values)
    # Add names via NumPy??
    # Leaving as R code for now.
    # names(values) <- .all_rows(data_frame)
    return values