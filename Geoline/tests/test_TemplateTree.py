from unittest import TestCase

from magby import Magby
from ..tests.testUtils import *

from ..Geoline.TemplateTree import *

url = 'https://magma.ucsf.edu'
token = 'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImFudG9uLm9nb3JvZG5pa292QHVjc2YuZWR1IiwibmFtZSI6IkFudG9uIE9nb3JvZG5pa292IiwicGVybSI6ImE6Y29wcm9qZWN0c190ZW1wbGF0ZSx4aGx0MTtlOmV4YW1wbGU7djpkc2NvbGFiLGlwaSxtdmlyMSIsImV4cCI6MTYxNDkxODQ1NX0.me0F4jm8Zwgndx5QTr9s0HqTzXI1SGMJvAn93EZ4OTrJBGGEib9gM9gf2ck1iQc4SSkNQOGTXLD-Tro9sYacaF8jgXKYotGcpWgQgWJHSqnYLcbMk9TtD3HRlgUQxRsldINI0bDDfr9tCdNoU4xz8Z52ydQ9r76rJY5Q-DomqbIO_4v7Ut0mKHljs7cVSwDWD5fT46tkKQqN6b-7gvGlmI9GceSHxfuTiW0juZNFEzJrcFK6AHZ2lWe7F_fSeoJdeHykEBx9kHXiTRGPDoQDO-N0zVdotEXtroROCWsPRmSc9O1G_yN7K7yFVaRwS-ghe43rAhJE4JrdKPWpHPa8sTUQp9Xf81IgEtkKuev3F5TLjP9gE20RxsLobnmO0ZypHsRciMRT4RdIBZ9xOT8Wmld97VQ-N9XMSJpkNnBoXOSwIo40YG5BrmO1EpTbiSbVwr3z9sKMd42LxSd9pf9s0xmiI7o7YX0yJir-QiBWYScUvx6qMR79XhfFubGvU7VDRvVVgxkIxgwDwW9B-VaauaazQvvMsSYmQNjGtyyOeNWAn_Uw88usxeN4qVOUvsS-1Niw2QZ_lNnSzt2QS5TVH-_PJrOfeEiIssjtzilAm57e09_hNYr_AptR9u1uCyF0VL_jiAAMke4wZBYFMx0D4SFJruD60WEAz0Cu2D-8e7Y'
project = 'example'
currDir = os.path.dirname(os.path.realpath(__file__))



class TestTemplateTree(TestCase):
    def setUp(self) -> None:
        mb = Magby.Magby(url, token)
        self.session = Session()
        self.vcr = prepCassette(self.session, os.path.join(currDir, 'fixtures/cassettes'))
        with self.vcr as vcr:
            vcr.use_cassette('TreeTemplate_setup')
            globalTemplate = mb.retrieve(project, 'all', [], 'all', dataType='json', session=self.session)
            self.templateTree = TemplateTree(globalTemplate)

    def test__ascendTree(self):
        rnaTree = self.templateTree._ascendTree(['rna_seq'])
        self.assertTrue(isinstance(rnaTree, list))
        self.assertEqual(rnaTree[2], 'subject')


    def test__commonRoot(self):
        primaryTree = self.templateTree._ascendTree(['rna_seq'])
        secondaryTree = self.templateTree._ascendTree(['flow'])
        commonRoot = self.templateTree._commonRoot(primaryTree, secondaryTree)
        self.assertTrue(isinstance(commonRoot, str))
        self.assertEqual(commonRoot, 'biospecimen')

    def test_traverseToModel(self):
        newPath = self.templateTree.traverseToModel('rna_seq', 'flow')
        self.assertTrue(isinstance(newPath, list))
        self.assertEqual(newPath, ['biospecimen', 'flow'])