from unittest import TestCase
from unittest.mock import patch
from pandas import DataFrame

from magby import Magby
from ..tests.testUtils import *

from ..Geoline.Geoline import *

url = 'https://magma.ucsf.edu'
token = 'token'
project = 'example'
currDir = os.path.dirname(os.path.realpath(__file__))

section = {
    'the': 'legolas:elf',
    'lord': 'gimli:dwarf',
    'of the': 'grishnakh:orc',
    'rings': 'samwise:hobbit'
}

class TestGeoline(TestCase):
    def setUp(self) -> None:
        self.geoline = Geoline(url, token, project)
        self.session = Session()
        self.vcr = prepCassette(self.session, os.path.join(currDir,'fixtures/cassettes'))

    def test__updater(self):
        with patch('builtins.input', side_effect=['y', 'y', 'grishnakh:mordor_orc', '0', 'gorbag:mordor_orc', 'STOP']):
            updatedSection = self.geoline._updater(section)
            self.assertTrue(isinstance(updatedSection, dict))
            self.assertEqual(updatedSection['of the'], 'gorbag:mordor_orc')

    def test__selectWorkflow(self):
        with patch('builtins.input', side_effect=['y', 'y', 'grishnakh:mordor_orc',
                                                  'y', '1', 'sample:treatment', 'n', # characteristics
                                                  'STOP']):
            samples = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
            self.assertTrue(isinstance(samples, dict))
            self.assertEqual(samples['organism'], 'grishnakh:mordor_orc')
            self.assertEqual(samples['processed data file'], 'rna_seq:gene_expression')


    def test__grouper(self):
        with patch('builtins.input', side_effect=['y', 'y', 'grishnakh:mordor_orc',
                                                  'y', '1', 'sample:treatment', 'n', # characteristics
                                                  'STOP']):
            samples = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
            modelGroups = self.geoline._grouper(samples)
            self.assertTrue(isinstance(samples, dict))
            self.assertEqual(modelGroups['grishnakh']['organism'], ['grishnakh', 'mordor_orc'])
            self.assertEqual(modelGroups['rna_seq']['processed data file'], ['rna_seq','gene_expression'])


    def test_constructMultiModelQuery(self):
        with patch('builtins.input', side_effect=['y', 'y', 'y',
                                                  'y', '1', 'rna_seq:fraction', 'n',  # characteristics
                                                  'STOP']):
            samples = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
            modelGroups = self.geoline._grouper(samples)
            query = self.geoline.constructMultiModelQuery(modelGroups, 'rna_seq')
            self.assertTrue(isinstance(query, dict))
            self.assertEqual(len(query), 12)
            self.assertEqual(query[2][1][0], 'EXAMPLE-HS12-WB1')







