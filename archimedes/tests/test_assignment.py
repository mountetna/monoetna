
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
            ret = str + ", yes'm!"
            return ret
        }

        @affirm = yesm("Looks like rain")
        '''
    )

    assert payload['affirm'] == "Looks like rain, yes'm!"
