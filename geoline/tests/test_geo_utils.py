import unittest
from unittest import TestCase
from unittest.mock import patch


from ..geoline.geo_utils import *

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

    def test_ask_attribute(self):
        with patch('builtins.input', side_effect=['sample:rna_seq_sample']):
            aw = ask_attribute(field='sample', magma_attr='sample:magma_sample')
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'sample:rna_seq_sample')
        with patch('builtins.input', side_effect=['y']):
            aw = ask_attribute(field='sample', magma_attr='sample:magma_sample')
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'y')
            self.assertNotEqual(aw, 'FFFF')

    def test_ask_characteristics(self):
        with patch('builtins.input', side_effect=['3', 'sample:rna_seq']):
            aw = ask_characteristics()
            self.assertTrue(isinstance(aw, tuple))
            self.assertEqual(aw, ('treatment', 'sample:rna_seq'))

    def test_add_another(self):
        with patch('builtins.input', side_effect=['STOP']):
            aw = add_another()
            self.assertTrue(isinstance(aw, str))
            self.assertEqual(aw, 'STOP')


    def test_characteristics(self):
        with patch('builtins.input', side_effect=['y', '1', 'model:cancer', 'n']):
            aw = characteristics(add_another, {})
            self.assertTrue(isinstance(aw, dict))
            self.assertEqual(aw, {'tissue': 'model:cancer'})



if __name__ == '__main__':
    unittest.main()