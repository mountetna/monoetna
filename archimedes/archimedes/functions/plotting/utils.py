import numpy as np

DISCRETE_KINDS = 'ObUS'
CONTINUOUS_KINDS = 'ifuc'

def _is_discrete(nda: np.ndarray):
    return nda.dtype.kind in DISCRETE_KINDS

def _is_continuous(nda: np.ndarray):
    return nda.dtype.kind in CONTINUOUS_KINDS

def _is_logical(nda: np.ndarray):
    return nda.dtype.kind == 'b'

def _is_integer(nda: np.ndarray):
    return nda.dtype.kind == 'i'

def _scale(y):
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

def _all_rows(df):
    return np.array(df.index)

def _which_rows(rows_use, data_frame):
    '''
    # converts a 'rows_use' given as string, logical, or numeric vector
    # into the df.index format expected internally by dittoViz functions.
    # Integer numeric vectors are treated as positional indexes, while string and other numeric types are treated as pd.DataFrame.index values
    # No matter the order given, the return from this function will be the targeted df.idex'es in the same order as they exist in the input data_frame.
    '''
    all_rows = _all_rows(data_frame)
    if rows_use==None:
        return all_rows
    
    r_use = np.array(rows_use)
    if is_logical(r_use):
        if (len(r_use)!=data_frame.shape[0]):
            raise Exception("'rows_use' length must equal the number of rows in 'data_frame' when given in logical form.")
        return all_rows[r_use]
    
    if is_integer(r_use):
        return(all_rows[r_use])
    
    # If here, string or numeric that should be of df.index
    return np.array([i for i in all_rows if i in r_use])
