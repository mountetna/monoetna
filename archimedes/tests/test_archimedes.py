import os
import sys

# setting working directory to the location of the test
curr_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_dir)
sys.path.append("../")

import tempfile
import pytest
from request_handler import resolve
import json
from errors import ArchimedesError
from pandas import Series

def test_reports_line_numbers():
    with pytest.raises(ArchimedesError) as e_info:
        payload = resolve('''@var1= 1
            @var2 = invalid syntax
            @var3 = 2''', '')

    assert e_info.value.info == 'Syntax error in line 2: No such variable invalid'


def test_supports_vectors():
    payload = resolve("@var1 = [ ant: 'a', bear: 'b', cat: 'c' ]")
    vector = payload['var1']
    assert type(vector) is Series
    assert vector.index.tolist() == ['ant', 'bear', 'cat']
    assert vector.values.tolist() == ['a', 'b', 'c']

def test_supports_indexing_into_vectors():
    payload = resolve(
     '''@list = [ ant: 'a', bear: 'b', cat: 'c' ]
        @int = @list[0]
        @let = @list['bear']'''
    )
    assert payload['int'] == 'a'
    assert payload['let'] == 'b'
 
def test_supports_nested_vectors():
    payload = resolve(
     '''@list1 = [ ant: 'a', bear: 'b', cat: 'c' ]
        @list2 = [ uno: 1, dos: 2, tres: 3 ]
        @list3 = [ list1: @list1, list2: @list2 ]'''
    )
    assert payload['list3']['list1'].tolist() == [ 'a', 'b', 'c' ]
    assert payload['list3']['list2'].tolist() == [ 1, 2, 3 ]

def test_supports_strings():
    payload = resolve("@var1 = 'a string'")
    assert payload['var1'] == 'a string'

def test_supports_numbers():
    payload = resolve('''
        @int = 12345
        @float = 1.2345
        @neg = -1.2345
        @sci = 1.23e-3
        ''')
    assert payload['int'] == 12345
    assert payload['float'] == 1.2345
    assert payload['neg'] == -1.2345
    assert payload['sci'] == 1.23e-3
