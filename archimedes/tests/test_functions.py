import os
import sys
import pytest

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append("../")

from errors import ArchimedesError
from helpers import resolve_json, labels, values

def test_defines_functions():
    payload = resolve_json('''
        yesm = (str) => {
            msg = str + ", yes'm!"
            return msg
        }

        @affirm = yesm("Looks like rain")
        ''')

    assert payload['affirm'] == "Looks like rain, yes'm!"

def test_defines_short_functions():
    payload = resolve_json('''
        yesm = (str) => str + ", yes'm!"

        @affirm = yesm("Looks like rain")
        ''')

    assert payload['affirm'] == "Looks like rain, yes'm!"

def test_allows_default_arguments():
    payload = resolve_json('''
        yesm = (str='Wonderful day') => {
            msg = str + ", yes'm!"
            return msg
        }

        @affirm = yesm()
        ''')

    assert payload['affirm'] == "Wonderful day, yes'm!"

def test_expects_argument_without_defaults():
    with pytest.raises(ArchimedesError) as e_info:
        payload = resolve_json('''
            yesm = (str) => {
                msg = str + ", yes'm!"
                return msg
            }

            @affirm = yesm()
            ''')

    assert e_info.value.info == 'Missing argument str with no default'

def test_allows_default_nil_arguments():
    payload = resolve_json('''
        yesm = (str=nil) => {
            msg = str ? str + ", yes'm!" : "Yes'm!"
            return msg
        }

        @affirm = yesm()
        ''')

    assert payload['affirm'] == "Yes'm!"

def test_maintains_separate_scopes():
    payload = resolve_json('''
        msg = 'blah'
        yesm = (str='Wonderful day') => {
            msg = str + ", yes'm!"
            return msg
        }

        @affirm = yesm()
        @message = msg
        ''')

    assert payload['affirm'] == "Wonderful day, yes'm!"
    assert payload['message'] == 'blah'

def test_reads_from_parent_scope():
    payload = resolve_json('''
        msg = 'blah'
        blah = () => msg

        @message = blah()
        ''')

    assert payload['message'] == 'blah'
