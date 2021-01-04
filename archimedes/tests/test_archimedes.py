import os
import sys
import pytest
from pandas import Series

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append("../")

from errors import ArchimedesError
from helpers import resolve_json, labels, values

def test_reports_line_numbers():
    with pytest.raises(ArchimedesError) as e_info:
        payload = resolve_json('''@var1= 1
            @var2 = invalid syntax
            @var3 = 2''')

    assert e_info.value.info == 'Syntax error in line 2: No such variable invalid'


def test_supports_vectors():
    payload = resolve_json("@var1 = [ ant: 'a', bear: 'b', cat: 'c' ]")
    vector = payload['var1']
    assert 'vector' in vector
    assert labels(vector) == ['ant', 'bear', 'cat']
    assert values(vector) == ['a', 'b', 'c']

def test_supports_indexing_into_vectors():
    payload = resolve_json(
     '''@list = [ ant: 'a', bear: 'b', cat: 'c' ]
        @int = @list[0]
        @let = @list['bear']'''
    )
    assert payload['int'] == 'a'
    assert payload['let'] == 'b'
 
def test_supports_nested_vectors():
    payload = resolve_json(
     '''@list1 = [ ant: 'a', bear: 'b', cat: 'c' ]
        @list2 = [ uno: 1, dos: 2, tres: 3 ]
        @list3 = [ list1: @list1, list2: @list2 ]'''
    )
    assert labels(payload['list3'])[0] == 'list1'
    assert values(values(payload['list3'])[0]) == [ 'a', 'b', 'c' ]

    assert labels(payload['list3'])[1] == 'list2'
    assert values(values(payload['list3'])[1]) == [ 1, 2, 3 ]

def test_supports_strings():
    payload = resolve_json("@var1 = 'a string'")
    assert payload['var1'] == 'a string'

def test_supports_numbers():
    payload = resolve_json('''
        @int = 12345
        @float = 1.2345
        @neg = -1.2345
        @sci = 1.23e-3
        ''')
    assert payload['int'] == 12345
    assert payload['float'] == 1.2345
    assert payload['neg'] == -1.2345
    assert payload['sci'] == 1.23e-3

def test_supports_local_vars():
    payload = resolve_json('''
        local = 2
        @ret = local + local ^ local
    '''
    )

    assert not 'local' in payload
    assert payload['ret'] == 6

def test_supports_access_operator():
    payload = resolve_json('''
        x = [ ant: 1, bear: 2, cat: 3 ]
        @y = x$bear
    ''')

    assert payload['y'] == 2
