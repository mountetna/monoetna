from unittest import TestCase
from unittest.mock import patch
from ..Geoline.geoUtils import *

testDict = {
    'a': 'lol',
    'b': {
        'x': 'kek',
        'y': 'gug'
    },
    'c': 'bab'
}

class Test(TestCase):
    def test_flatten(self):
        flat = flatten(testDict, '_')
        self.assertTrue(isinstance(flat, dict))
        self.assertEqual(list(flat.keys())[1], 'b_x')
        self.assertEqual(list(flat.values())[2], 'gug')

    def test_askAttribute(self):
        with patch('builtins.input', side_effect=['rna_seq_sample']):
            aw = askAttribute(field='sample', magmaAttr='magma_sample')
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'rna_seq_sample')
        with patch('builtins.input', side_effect=['Y']):
            aw = askAttribute(field='sample', magmaAttr='magma_sample')
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'Y')
            self.assertNotEqual(aw, 'FFFF')

    def test_askCharacteristics(self):
        with patch('builtins.input', side_effect=['3', 'rna_seq']):
            aw = askCharacteristics()
            self.assertTrue(isinstance(aw, tuple))
            self.assertEqual(aw, ('treatment', 'rna_seq'))

    def test_addAnother(self):
        with patch('builtins.input', side_effect=['STOP']):
            aw = addAnother()
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'STOP')
