import pandas as pd
import numpy as np
from typing import Any, Callable

def _can_base(x: Any, converter: Callable[[Any], Any], na_ok=True):
    x_series = pd.Series(x)
    coerced = x_series.apply(converter)
    result = ~coerced.isna()
    if na_ok:
        result = result | x_series.isna()
    if len(result)>1:
        return result.tolist()
    return result.tolist()[0]

def _can_numeric(x: Any, na_ok=True):
    return _can_base(x, lambda x: pd.to_numeric(x, errors='coerce'), na_ok)

def _to_logical(val: Any):
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        val_lower = val.strip().lower()
        if val_lower in ['true', 't', 'y', 'yes']:
            return True
        elif val_lower in ['false', 'f', 'n', 'no']:
            return False
        else:
            return np.nan
    if pd.isna(val):
        return np.nan
    if _can_numeric(val):
        if val in [0,1]:
            return bool(val)
        else:
            return np.nan
    raise Exception(f'Unaccounted for data-type in _to_logical, {val}, of type {type(val)}')

def _can_logical(x: Any, na_ok=True):
    return _can_base(x, _to_logical, na_ok)

def _to_str(val: Any):
    if isinstance(val, str):
        return val
    if pd.isna(val):
        return np.nan
    try:
        return str(val)
    except:
        return np.nan

def _can_str(x: Any, na_ok=True):
    return _can_base(x, lambda x: _to_str(x), na_ok)