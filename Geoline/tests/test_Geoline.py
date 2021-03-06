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


    def test__constructMultiModelQuery(self):
        with patch('builtins.input', side_effect=['y', 'y', 'y',
                                                  'y', '1', 'rna_seq:fraction', 'n',  # characteristics
                                                  'STOP']):
            samples = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
            modelGroups = self.geoline._grouper(samples)
            with self.vcr as vcr:
                vcr.use_cassette('Geoline__constructMultiModelQuery')
                query = self.geoline._constructMultiModelQuery(modelGroups, 'rna_seq', session=self.session)
            self.assertTrue(isinstance(query, list))
            self.assertEqual(len(query), 3)
            self.assertEqual(query[2][4], 'gene_expression')


    def test__walkAnswer(self):
        with patch('builtins.input', side_effect=['y', 'flow:stain', 'subject:group',
                                                  'y', '1', 'rna_seq:fraction', 'n',  # characteristics
                                                  '', '', '']):
            samples = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
            modelGroups = self.geoline._grouper(samples)
            with self.vcr as vcr:
                vcr.use_cassette('Geoline__walkAnswer')
                query = self.geoline._constructMultiModelQuery(modelGroups, 'rna_seq', session=self.session)
                mb = Magby(url, token)
                reply = mb.query(self.geoline._projectName, query, session=self.session)
                aw = self.geoline._walkAnswer(reply['answer'][0], reply['format'])
            self.assertTrue(isinstance(aw, dict))
            self.assertEqual(aw['rna_seq:tube_name'], 'EXAMPLE-HS10-WB1-RSQ1')
            self.assertEqual(len(aw), 5)

    def test__walkAnswerWithOneToMany(self):
        awElem = ['SAMPLE1-RSQ1',
                  ['Lembas Bread',
                   [['Bread Sort 1', None],
                    ['Bread Sort 3', 'Honey'],
                    ['Bread Sort 5', None]]]]
        awFormat = ['proj::elven_bread#shape',
                    ['proj::elven_bread#texture',
                     ['proj::elven_bread#pool_name', 'proj::bread#panel']]]

        aw = self.geoline._walkAnswer(awElem, awFormat)
        self.assertTrue(isinstance(aw, dict))
        self.assertEqual(aw['elven_bread:pool_name'], 'Bread Sort 1; Bread Sort 3; Bread Sort 5')


    def test__reduceOneToMany(self):
        otm = [['Name1', None],
               ['Name2', 'Bread Sample v1'],
               ['Name3', None]]
        reduced = self.geoline._reduceOneToMany(otm)
        self.assertTrue(isinstance(reduced[0], str))
        self.assertEqual(reduced[1], 'None; Bread Sample v1; None')







