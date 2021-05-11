import unittest
from unittest import TestCase
from ..magby.Magby import *
from .testUtils import *


url = 'https://magma.ucsf.edu'
token = 'token'
currDir = os.path.dirname(os.path.realpath(__file__))


class TestMagby(TestCase):
    def setUp(self) -> None:
        self.session = Session()
        self.vcr = prepCassette(self.session, os.path.join(currDir,'fixtures/cassettes'))
        self.magby = Magby(url, token)

    def test_url(self):
        self.assertEqual(self.magby.url, url)
        self.magby.url = 'changed_url'
        self.assertNotEqual(self.magby.url, url)
        self.assertEqual(self.magby.url, 'changed_url')

    def test_token(self):
        self.assertEqual(self.magby.token, token)
        self.magby.token = 'changed_token'
        self.assertNotEqual(self.magby.token, token)
        self.assertEqual(self.magby.token, 'changed_token')

    def test__constructPayload(self):
        payload = self.magby._constructPayload(projectName='example',
                                               modelName='subject',
                                               recordNames=["EXAMPLE-HS1"],
                                               attributeNames=["group"],
                                               formatBy='json')
        self.assertTrue(isinstance(payload, dict))


    def test_getProjects(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_getProjects')
            projects = self.magby.getProjects(session=self.session)
        self.assertEqual(projects[2]['project_name'], 'ipi')
        self.assertTrue(isinstance(projects, list))
        self.assertTrue(len(projects) > 0)

    def test_retrieve(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_retrieve')
            out = self.magby.retrieve(projectName="example", modelName='subject',
                                      recordNames='all', attributeNames="group", dataType='meta',
                                      session=self.session)
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (1618,2))

    def test_retrieveJSON(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_retrieveJSON')
            out = self.magby.retrieve(projectName="example", modelName='subject',
                                      recordNames='all', attributeNames="group", dataType='json',
                                      session=self.session)
        self.assertTrue(isinstance(out, dict))
        self.assertEqual(out['models']['subject']['documents']["EXAMPLE-HS1"]['group'], 'g1')
        self.assertEqual(len(out['models']['subject']['documents']), 12)


    def test_retrieveMatrix(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_retrieveMatrix')
            out = self.magby.retrieve(projectName="example", modelName='rna_seq',
                                      recordNames='all', attributeNames=["gene_tpm"], dataType='mtx',
                                      session=self.session)
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (40,12))
        self.assertRaises(MagmaError, self.magby.retrieve, projectName="example", modelName='rna_seq',
                          recordNames='all', attributeNames=["gene_counts", "gene_tpm"],
                          dataType='mtx',
                          session=self.session)


    def test_query(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_query')
            out = self.magby.query(projectName='example', queryTerms=['rna_seq', '::all', 'biospecimen', '::identifier'],
                                   session=self.session)
        self.assertTrue(isinstance(out, dict))
        self.assertEqual(out['answer'][0][1], 'EXAMPLE-HS10-WB1')
        self.assertEqual(len(out['answer']), 12)


if __name__ == '__main__':
    unittest.main()
