import os
import sys
import pytest

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append("../")

from errors import ArchimedesError
from helpers import resolve_json, labels, values

def test_allows_true_false_and_nil_keywords():
    payload = resolve_json(
     '''@t = true
        @f = false
        @n = nil'''
    )
    assert payload['t'] == True
    assert payload['f'] == False
    assert payload['n'] == None

def test_it_compares_two_numbers():
    payload = resolve_json('@vec = 10 > 2')
    assert payload['vec'] == True

def test_compares_two_vectors():
    payload = resolve_json(
     '''@vec1 = [ 1, 2, 3 ]
      @vec2 = [ 1, 4, 3 ]
      @same = @vec1 == @vec2'''
    )
    assert values(payload['same']) == [ True, False, True ]

def test_compares_a_vector_and_a_number():
    payload = resolve_json(
     '''@vec1 = [ 1, 2, 3 ]
      @gte = @vec1 >= 2'''
    )
    assert values(payload['gte']) == [ False, True, True ]

def test_compares_to_nil():
    payload = resolve_json(
     '''@vec1 = [ 1, nil, 3 ]
      @n1 = @vec1 == nil
      @n2 = @vec1 != nil'''
    )
    assert values(payload['n1']) == [ False, True, False ]
    assert values(payload['n2']) == [ True, False, True ]

def test_supports_boolean_comparisons():
    payload = resolve_json(
        '''@or1 = true || false
        @or2 = false || false
        @or3 = false || true
        @or4 = true || true
        @and1 = true && false
        @and2 = false && false
        @and3 = false && true
        @and4 = true && true'''
    )
    assert payload['or1'] == True
    assert payload['or2'] == False
    assert payload['or3'] == True
    assert payload['or4'] == True
    assert payload['and1'] == False
    assert payload['and2'] == False
    assert payload['and3'] == False
    assert payload['and4'] == True

def test_supports_boolean_vector_comparisons():
    payload = resolve_json(
      '''@v1 = [ true, false, false, true ]
       @v2 = [ false, false, true, true ]
       @or = @v1 || @v2
       @and = @v1 && @v2'''
    )
    assert values(payload['or']) == [ True, False, True, True ]
    assert values(payload['and']) == [ False, False, False, True ]

def test_supports_math_operations():
    payload = resolve_json(
      '''
      @mul = 4 * 4
      @div = 4 / 4
      @add = 4 + 4
      @exp = 4 ^ 4
      @sub = 4 - 4
      @ternary = 4 == 4 ? 4 : -4
      '''
    )
    assert payload['mul'] == 16
    assert payload['div'] == 1
    assert payload['add'] == 8
    assert payload['exp'] == 256
    assert payload['sub'] == 0
    assert payload['ternary'] == 4

def test_obeys_precedence():
    payload = resolve_json(
      '''
      @mul_plus = 4 * 4 + 2
      @div_sub = 4 / 4 - 2
      @add_exp = 4 + 4 ^ 2
      '''
    )
    assert payload['mul_plus'] == 18
    assert payload['div_sub'] == -1
    assert payload['add_exp'] == 20

def test_supports_infix_notation():
    payload = resolve_json( '@calc = (4 + 4) * 4 + 4')

    assert payload['calc'] == 36

def test_supports_comparisons():
    payload = resolve_json(
     '''@gt = 4 > 0
      @gte = 4 >= 9
      @lt = 4 < 10
      @lte = 4 <= 15
      @eq = 4 == 15'''
    )

    assert payload['gt'] == True
    assert payload['gte'] == False
    assert payload['lt'] == True
    assert payload['lte'] == True
    assert payload['eq'] == False

def test_binds_a_vector_of_column_vectors_into_a_matrix():
    payload = resolve_json('''@x = :: [ a: [ i: 1, j: 2 ], b: [ 3, 4 ] ]''')

    assert payload['x'] == {
        'matrix': {
            'col_names': ['a', 'b'],
            'row_names': ['i', 'j'],
            'rows': [[1, 3], [2, 4]]
        }
    }

def test_binds_to_an_existing_matrix():
    payload = resolve_json('''x = :: [ a: [ i: 1, j: 2 ], b: [ 3, 4 ] ]
        @y = x :: [ c: [ 5, 6 ] ]
        ''')

    assert payload['y'] == {
        'matrix': {
            'col_names': ['a', 'b', 'c'],
            'row_names': ['i', 'j'],
            'rows': [[1, 3, 5], [2, 4, 6]]
        }
    }
