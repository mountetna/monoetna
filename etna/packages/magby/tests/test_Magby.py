from unittest import TestCase
import configparser
from requests import Session
from ..magby.Magby import *
from ..tests.testUtils import *

conf = configparser.ConfigParser()
conf.read('magby/tests/proxyConfig.ini')    # Config is in .gitignore



class TestMagby(TestCase):
    def setUp(self) -> None:
        self.proxy = json.loads(conf['DEFAULT'].get('proxy'))
        self.session = Session()
        self.session.proxies.update(self.proxy)
        self.session.verify = False
        self.vcr = prepCassette(self.session, './magby/tests/fixtures/cassettes')
        self.magby = Magby(conf['DEFAULT'].get('url'),
                           conf['DEFAULT'].get('token'))

    def test_url(self):
        self.assertEqual(self.magby.url, conf['DEFAULT'].get('url'))
        self.magby.url = 'changed_url'
        self.assertNotEqual(self.magby.url, conf['DEFAULT'].get('url'))
        self.assertEqual(self.magby.url, 'changed_url')

    def test_token(self):
        self.assertEqual(self.magby.token, conf['DEFAULT'].get('token'))
        self.magby.token = 'changed_token'
        self.assertNotEqual(self.magby.token, conf['DEFAULT'].get('token'))
        self.assertEqual(self.magby.token, 'changed_token')

    def test__constructPayload(self):
        payload = self.magby._constructPayload(projectName='ipi',
                                               modelName='sample',
                                               recordNames=["IPIADR001.T1"],
                                               attributeNames=["patient"],
                                               format='json')
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
            out = self.magby.retrieve(projectName="ipi", modelName='sample',
                                      recordNames='all', attributeNames=["patient"], dataType='meta',
                                      session=self.session)
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (771,2))

    def test_retrieveJSON(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_retrieveJSON')
            out = self.magby.retrieve(projectName="ipi", modelName='sample',
                                      recordNames='all', attributeNames=["patient"], dataType='json',
                                      session=self.session)
        self.assertTrue(isinstance(out, dict))
        self.assertEqual(out['models']['sample']['documents']['IPIADR001.T1']['patient'], 'IPIADR001')
        self.assertEqual(len(out['models']['sample']['documents']), 771)


    def test_retrieveMatrix(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_retrieveMatrix')
            out = self.magby.retrieve(projectName="ipi", modelName='rna_seq',
                                      recordNames=['IPIBLAD005.N1.rna.cd45neg'], attributeNames=["gene_counts"], dataType='mtx',
                                      session=self.session)
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (58051,1))
        self.assertRaises(MagmaError, self.magby.retrieve, projectName="ipi", modelName='rna_seq',
                          recordNames=['IPIBLAD005.N1.rna.cd45neg'], attributeNames=["gene_counts", "gene_tpm"],
                          dataType='mtx',
                          session=self.session)


    def test_query(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magby_query')
            out = self.magby.query(projectName='ipi', queryTerms=['sample', '::all', 'patient', '::identifier'],
                                   session=self.session)
        self.assertTrue(isinstance(out, dict))
        self.assertEqual(out['answer'][0][1], 'IPIADR001')
        self.assertEqual(len(out['answer']), 771)



