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
    '''
    Centers and scales a vector of numeric values to have mean=0 and stdev=1
    'y' and numeric vector of at least length > 1.
    '''
    x = y.copy()
    x -= x.mean()
    x /= x.std()
    return x

def _leave_default_or_null(target, default, none_if = False, default_when = "make"):
    '''
    A system for filling in input defaults whenever a given 'target' input is left with the dafaut value given within the function definition.
    
    Necessary because one cannot have param2 have a default value of param1 in python function definitions.
    Useful for many use cases.
    
    'target' the input value
    'default' truly intended default value
    'default_when' the default value used within the funciton definition
    'none_if' a bool which, when True, will cause the function to output 'None' instead of 'default'
    
    Value: Either 'target' itself, 'default', or literally 'None'.
    
    Example:
    # In a function definition where 'xlab' is intended to become whatever the user gave to 'x_by',
    #  unless the user explicitly provides something different to 'xlab'...
    # Function definition should have 'xlab = "make"'
    # Then, inside the function code:
    xlab = _leave_default_or_null(x_lab, x_by)
    '''
    if target == default_when:
        if none_if:
            return None
        else:
            return default
    return target

def _all_rows(df):
    '''
    Returns the standard identifiers of 'df' row identifiers that Viz fundtions will rely on.
    Centralized function abstraction makes it easier to update row indexing if needed.
    '''
    return df.index

def _which_rows(rows_use, data_frame):
    '''
    The workhorse of the flexible rows selection system.
    Converts a 'rows_use' given as string list, logical list, or numeric list into the 'data_frame.index' format used directly by Viz functions.
    
    'rows_use' A string list, logical list, or numeric list
    'data_frame' a pandas DataFrame
    
    Value: (subset of) data_frame.index targeted by 'rows_use'.
    Note: Integer numeric vectors are treated as positional indexes, while string and other numeric types are expected to be pd.DataFrame.index values.
    Note: The return will always be in the same order as indexes exist in data_frame.index, regardless of the order of 'rows_use' values.
    '''
    all_rows = _all_rows(data_frame)
    if rows_use==None:
        return all_rows
    
    r_use = np.array(rows_use)
    if _is_logical(r_use):
        if (len(r_use)!=data_frame.shape[0]):
            raise Exception("'rows_use' length must equal the number of rows in 'data_frame' when given in logical form.")
        return all_rows[r_use]
    
    if _is_integer(r_use):
        return(all_rows[r_use])
    
    # If here, string or numeric that should be of df.index
    return np.array([i for i in all_rows if i in r_use])
