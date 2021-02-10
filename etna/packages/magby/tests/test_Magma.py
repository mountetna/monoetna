from unittest import TestCase
import vcr
import configparser
from ..magby.Magma import *



conf = configparser.ConfigParser()
conf.read('magby/tests/proxyConfig.ini')    # Config is in .gitignore


class TestMagma(TestCase):
    @vcr.use_cassette('tests/fixtures/vcr_cassettes/test_Magma.yaml')
    def setUp(self) -> None:
        proxy = json.loads(conf['DEFAULT'].get('proxy'))
        self.magma = Magma(conf['DEFAULT'].get('url'),
                           conf['DEFAULT'].get('token'),
                           'retrieve')
        self.magma._session.proxies.update(proxy)

    def test_getResponseContent(self):
        payload = {
            "project_name": "ipi",
            "model_name": "sample",
            "record_names": ["IPIADR001.T1"],
            "attribute_names": ["patient"]
        }
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
        magmaResponse, magmaResponseHeaders = self.magma.magmaCall(payload)
        self.assertTrue(isinstance(magmaResponse, dict))
        self.assertTrue(isinstance(magmaResponseHeaders, dict))
        self.assertEqual(magmaResponse['models']['sample']['documents']["IPIADR001.T1"]['patient'], 'IPIADR001')
        self.assertEqual(magmaResponseHeaders, 'headers_stuff') # TODO get actual data






    
