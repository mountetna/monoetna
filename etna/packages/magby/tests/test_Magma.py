from unittest import TestCase
import configparser
from ..magby.Magma import *
from ..tests.testUtils import *

conf = configparser.ConfigParser()
conf.read('magby/tests/proxyConfig.ini')    # Config is in .gitignore



class TestMagma(TestCase):

    def setUp(self) -> None:
        proxy = json.loads(conf['DEFAULT'].get('proxy'))
        self.magma = Magma(conf['DEFAULT'].get('url'),
                           conf['DEFAULT'].get('token'),
                           'retrieve')
        self.magma._session.proxies.update(proxy)
        self.vcr = prepCassette(self.magma._session)


    def test_getResponseContent(self):
        payload = {
            "project_name": "ipi",
            "model_name": "sample",
            "record_names": ["IPIADR001.T1"],
            "attribute_names": ["patient"]
        }
        with self.vcr as vcr:
            vcr.use_cassette('Mamgma_getResponseContent')
            response = self.magma._session.post(self.magma._url, data=json.dumps(payload), headers=self.magma._headers, verify=False)
            responseJSON = self.magma.getResponseContent(response)
        self.assertTrue(isinstance(responseJSON, dict))
        self.assertEqual(len(responseJSON['models']['sample']['documents']["IPIADR001.T1"].keys()), 2)
        self.assertEqual(responseJSON['models']['sample']['documents']["IPIADR001.T1"]['patient'], 'IPIADR001')


    def test_magmaCall(self):
        payload = {
            "project_name": "ipi",
            "model_name": "sample",
            "record_names": ["IPIADR001.T1"],
            "attribute_names": ["patient"]
        }
        with self.vcr as vcr:
            vcr.use_cassette('Mamgma_magmaCall')
            magmaResponse, magmaResponseHeaders = self.magma.magmaCall(payload, verify=False)
        self.assertTrue(isinstance(magmaResponse, dict))
        self.assertEqual(magmaResponse['models']['sample']['documents']["IPIADR001.T1"]['patient'], 'IPIADR001')
        self.assertEqual(magmaResponseHeaders['Content-Length'], '9018')



