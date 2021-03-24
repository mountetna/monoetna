import unittest
from unittest import TestCase
from magby import Magby

from testUtils import *
from geoline.TemplateTree import *

url = 'https://magma.ucsf.edu'
token = 'token'
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
        self.assertEqual(newPath, ['biospecimen', 'flow', '::all'])


if __name__ == '__main__':
    unittest.main()