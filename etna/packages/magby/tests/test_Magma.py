import unittest
from unittest import TestCase
from ..magby.Magma import *
from .testUtils import *

url = 'https://magma.ucsf.edu'
token = 'token'
currDir = os.path.dirname(os.path.realpath(__file__))

payload = {
    "project_name": "example",
    "model_name": "subject",
    "record_names": ["EXAMPLE-HS1"],
    "attribute_names": ["group"]
}

class TestMagma(TestCase):

    def setUp(self) -> None:
        self.magma = Magma(url, token, 'retrieve')
        self.vcr = prepCassette(self.magma._session, os.path.join(currDir,'fixtures/cassettes'))


    def test_getResponseContent(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magma_getResponseContent')
            response = self.magma._session.post(self.magma._url, data=json.dumps(payload), headers=self.magma._headers, verify=False)
            responseJSON = self.magma.getResponseContent(response)
        self.assertTrue(isinstance(responseJSON, dict))
        self.assertEqual(len(responseJSON['models']['subject']['documents']["EXAMPLE-HS1"].keys()), 2)
        self.assertEqual(responseJSON['models']['subject']['documents']["EXAMPLE-HS1"]['group'], 'g1')


    def test_magmaCall(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magma_magmaCall')
            magmaResponse, magmaResponseHeaders = self.magma.magmaCall(payload)
        self.assertTrue(isinstance(magmaResponse, dict))
        self.assertEqual(magmaResponse['models']['subject']['documents']["EXAMPLE-HS1"]['group'], 'g1')
        self.assertEqual(magmaResponseHeaders['Content-Length'], '1376')

    def test__janusCall(self):
        with self.vcr as vcr:
            vcr.use_cassette('Magma__janusCall')
            janusProjects = self.magma._janusProjectsCall('/'.join([url.replace('magma', 'janus'),
                                                                    'projects']))
        self.assertEqual(janusProjects['projects'][2]['project_name'], 'ipi')
        self.assertTrue(isinstance(janusProjects, dict))
        self.assertTrue(len(janusProjects) > 0)


if __name__ == '__main__':
    unittest.main()



