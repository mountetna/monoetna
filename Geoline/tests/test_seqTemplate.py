from unittest import TestCase
from unittest.mock import patch

from ..Geoline.seqTemplate import characteristics
from ..Geoline.geoUtils import *

class Test(TestCase):
    def test_characteristics(self):
        with patch('builtins.input', side_effect=['y', '1', 'cancer', 'n']):
            aw = characteristics(addAnother, {})
            self.assertTrue(isinstance(aw, dict))
            self.assertEqual(aw, {'tissue': 'cancer'})

