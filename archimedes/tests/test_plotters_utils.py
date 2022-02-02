from archimedes.functions.plotting.utils import _all_rows, _leave_default_or_none, _which_rows, _is_continuous, _is_discrete, _is_integer, _is_logical
from tests.fixtures.example_df import example_df
import numpy as np

# Columns of example_df used
cont1 = '0'
disc1 = 'leiden'

str_arr = np.array(["A","B","C"])
bool_arr = np.array([True, False, False])
int_arr = np.array([1,2,3])
float_arr = np.array([1.5,2,3])

def test_is_discrete():
    assert _is_discrete(str_arr)
    assert _is_discrete(bool_arr)
    assert not _is_discrete(int_arr)
    assert not _is_discrete(float_arr)

def test_is_continuous():
    assert not _is_continuous(str_arr)
    assert not _is_continuous(bool_arr)
    assert _is_continuous(int_arr)
    assert _is_continuous(float_arr)

def test_is_logical():
    assert not _is_logical(str_arr)
    assert _is_logical(bool_arr)
    assert not _is_logical(int_arr)
    assert not _is_logical(float_arr)

def _is_integer():
    assert not _is_integer(str_arr)
    assert not _is_integer(bool_arr)
    assert _is_integer(int_arr)
    assert not _is_integer(float_arr)

def test_all_rows():
    df = example_df().copy()
    all_rows = _all_rows(df)
    
    assert df.loc[all_rows].shape == df.shape

def test_which_rows():
    df = example_df().copy()
    all_rows = _all_rows(df)
    
    # None = _all_rows
    which_rows_none = _which_rows(None, df)
    assert df.loc[which_rows_none].shape == df.shape
    assert all(which_rows_none == all_rows)
    
    # bool iterable = True indexes
    logic = [ df[cont1][i]>5 for i in range(df.shape[0]) ]
    
    which_rows_logic = _which_rows(logic, df)
    assert df.loc[which_rows_logic].shape[0] == sum(logic)
    assert all(i in all_rows for i in which_rows_logic)
    
    # int iterable = location based indexes
    ints = [0,1,2,3,10]
    
    which_rows_ints = _which_rows(ints, df)
    should_match = all_rows[ints]
    assert df.loc[which_rows_ints].shape[0] == len(ints)
    assert all(which_rows_ints == should_match)
    
    # other interable = expected as values df.index in df.index order
    indexes = ['0','1','2','3','10']
    reordered_indexes = indexes[::-1]
    
    which_rows_indexes = _which_rows(indexes, df)
    assert df.loc[which_rows_indexes].shape[0] == len(indexes)
    assert all(which_rows_indexes == _which_rows(reordered_indexes, df))

def test_leave_default_or_null():
    assert _leave_default_or_none(None, "default")==None
    assert _leave_default_or_none("make", "default")=="default"
    assert _leave_default_or_none("make", "default", none_if=True)==None
    assert _leave_default_or_none("make", "default", default_when="not_me")=="make"

