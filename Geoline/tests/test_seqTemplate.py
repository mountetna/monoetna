from unittest import TestCase
import os
from ..Geoline.seqTemplate import characteristics
from ..Geoline.geoUtils import *

class Test(TestCase):
    def test_characteristics(self):
        print(os.getcwd())
        a = characteristics(addAnother, {})
        pass
