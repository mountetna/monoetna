import os
import sys
import pytest

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append("../")

from errors import ArchimedesError
from helpers import resolve_json, labels, values

def test_handles_builtins():
    payload = resolve_json('''
        vec = [ 1, 2, 3, 4 ]

        @length = length(vec)
        ''')

    assert payload['length'] == 4

def test_expects_arguments_for_builtins():
    with pytest.raises(ArchimedesError) as e_info:
        payload = resolve_json('''
            @length = length()
            ''')

    assert e_info.value.info == "length() missing 1 required positional argument: 'vec'"


def test_max():
    payload = resolve_json('''
        vec = [ 1, 2, 3, 4 ]

        @max = max(vec)
        ''')

    assert payload['max'] == 4

def test_min():
    payload = resolve_json('''
        vec = [ 1, 2, 3, 4 ]

        @min = min(vec)
        ''')

    assert payload['min'] == 1

def test_unique():
    payload = resolve_json('''
        vec = [ 1, 1, 3, 'a', 'a' ]

        @unique = unique(vec)
        ''')

    assert values(payload['unique']) == [ 1, 3, 'a' ]

def test_join():
    payload = resolve_json('''
        vec = [ 1, 1, 3, 'a', 'a' ]

        @join = join(vec)
        ''')

    assert payload['join'] == '113aa'

def test_rep():
    payload = resolve_json('''
        vec = [ 1 ]

        @rep = rep(vec, 3)
        ''')

    assert values(payload['rep']) == [ 1, 1, 1 ]

def test_seq():
    payload = resolve_json('''
        @seq = seq(1, 5, 2)
        ''')

    assert values(payload['seq']) == [ 1, 3, 5 ]

def test_compact():
    payload = resolve_json('''
        @vec = compact( [1, 2, nil, 3, nil ] )
        ''')

    assert values(payload['vec']) == [ 1, 2, 3 ]

def test_labels():
    payload = resolve_json('''
        @labels = labels([ a: 1, b: 2, c: 3 ])
        ''')

    assert values(payload['labels']) == [ 'a', 'b', 'c' ]

def test_set_labels():
    payload = resolve_json('''
        @vec = set_labels([ 1, 2, 3 ], [ 'a', 'b', 'c' ])
        ''')

    assert values(payload['vec']) == [ 1, 2, 3 ]
    assert labels(payload['vec']) == [ 'a', 'b', 'c' ]

def test_filter():
    payload = resolve_json('''
        @vec = filter([ 1, 2, 3, 4, 5 ], (x) => x >= 3)
        ''')

    assert values(payload['vec']) == [ 3, 4, 5 ]

def test_any():
    payload = resolve_json('''
        @any = any([ 1, 2, 3, 4, 5 ], (x) => x >= 3)
        ''')

    assert payload['any'] == True

def test_all():
    payload = resolve_json('''
        @all = all([ 1, 2, 3, 4, 5 ], (x) => x >= 3)
        ''')

    assert payload['all'] == False

def test_group():
    payload = resolve_json('''
        roster = [
            [ day: 'Saturday', name: 'Thales' ],
            [ day: 'Sunday', name: 'Thales' ],
            [ day: 'Monday', name: 'Eudoxus' ],
            [ day: 'Tuesday', name: 'Thales' ],
            [ day: 'Wednesday', name: 'Eudoxus' ],
            [ day: 'Thursday', name: 'Thales' ],
            [ day: 'Friday', name: 'Philolaus' ],
            [ day: 'Saturday', name: 'Eudoxus' ]
        ]
        @duties = group(roster, (assignment) => assignment['name'])
        ''')

    assert labels(payload['duties']) == [ 'Thales', 'Eudoxus', 'Philolaus' ]

def test_rowsums():
    payload = resolve_json('''
        mat = ::[
            [ 1, 1, 1 ],
            [ 2, 2, 2 ],
            [ 3, 3, 3 ]
        ]
        @sum = rowsums(mat)
        ''')

    assert values(payload['sum']) == [ 3, 6, 9 ]

def test_colsums():
    payload = resolve_json('''
        mat = ::[
            [ 1, 1, 1 ],
            [ 2, 2, 2 ],
            [ 3, 3, 3 ]
        ]
        @sum = colsums(mat)
        ''')

    assert values(payload['sum']) == [ 6, 6, 6 ]

def test_rownames():
    payload = resolve_json('''
        mat = ::[
            a: [ i: 1, j: 1, k: 1 ],
            b: [ 2, 2, 2 ],
            c: [ 3, 3, 3 ]
        ]
        @names = rownames(mat)
        ''')

    assert values(payload['names']) == [ 'i', 'j', 'k' ]

def test_colnames():
    payload = resolve_json('''
        mat = ::[
            a: [ i: 1, j: 1, k: 1 ],
            b: [ 2, 2, 2 ],
            c: [ 3, 3, 3 ]
        ]
        @names = colnames(mat)
        ''')

    assert values(payload['names']) == [ 'a', 'b', 'c' ]

def test_numrows():
    payload = resolve_json('''
        mat = ::[
            a: [ i: 1, j: 1, k: 1 ],
            c: [ 3, 3, 3 ]
        ]
        @num = numrows(mat)
        ''')

    assert payload['num'] == 3

def test_numcols():
    payload = resolve_json('''
        mat = ::[
            a: [ i: 1, j: 1, k: 1 ],
            c: [ 3, 3, 3 ]
        ]
        @num = numcols(mat)
        ''')

    assert payload['num'] == 2
