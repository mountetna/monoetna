import unittest
from unittest import TestCase
from magby import Magby

from .test_utils import *
from ..geoline.TemplateTree import TemplateTree

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
            global_template = mb.retrieve(project, 'all', [], 'all', dataType='json', session=self.session)
            self.template_tree = TemplateTree(global_template)

    def test__ascend_tree(self):
        rna_tree = self.template_tree._ascend_tree(['rna_seq'])
        self.assertTrue(isinstance(rna_tree, list))
        self.assertEqual(rna_tree[2], 'subject')

    def test__common_root(self):
        primary_tree = self.template_tree._ascend_tree(['rna_seq'])
        secondary_tree = self.template_tree._ascend_tree(['flow'])
        common_root = self.template_tree._common_root(primary_tree, secondary_tree)
        self.assertTrue(isinstance(common_root, str))
        self.assertEqual(common_root, 'biospecimen')

    def test_traverse_to_model(self):
        new_path = self.template_tree.traverse_to_model('rna_seq', 'flow')
        self.assertTrue(isinstance(new_path, list))
        self.assertEqual(new_path, ['biospecimen', 'flow', '::all'])


unittest.main()
