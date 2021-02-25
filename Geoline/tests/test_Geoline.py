from unittest import TestCase
from unittest.mock import patch

from magby import Magby
from ..tests.testUtils import *

from ..Geoline.Geoline import *

url = 'https://magma.ucsf.edu'
token = 'token'
currDir = os.path.dirname(os.path.realpath(__file__))

section = {
    'the': 'legolas:elf',
    'lord': 'gimli:dwarf',
    'of the': 'grishnakh:orc',
    'rings': 'sam:hobbit'
}

class TestGeoline(TestCase):
    def setUp(self) -> None:
        self.geoline = Geoline(url, token)
        self.session = Session()
        self.vcr = prepCassette(self.session, os.path.join(currDir,'fixtures/cassettes'))

    def test__updater(self):
        with patch('builtins.input', side_effect=['y', 'y', 'grishnakh:mordor_orc', '0', 'gorbag:mordor_orc', 'STOP']):
            updatedSection = self.geoline._updater(section)
            self.assertTrue(isinstance(updatedSection, dict))
            self.assertEqual(updatedSection['of the'], 'gorbag:mordor_orc')

    def test__selectWorkflow(self):
        #with patch('builtins.input', side_effect=['y', 'y', 'grishnakh:mordor_orc',
        #                                          'y', '1', 'cancer', 'n', # characteristics
        #                                          'STOP']):
        sampleSection = self.geoline._selectWorkflow(samplesSection, 'rna_seq')
        pass




