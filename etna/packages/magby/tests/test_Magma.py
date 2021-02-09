from unittest import TestCase
from ..magby.Magma import *


class TestMagma(TestCase):
    def setUp(self) -> None:
        self.magma = Magma('some_url', 'my_token')
        self.json_string = '{"projects":[{"project_name":"mvir1","project_name_full":"COMET","project_description":null},' \
                           '{"project_name":"dscolab","project_name_full":"Data Science CoLab","project_description":null}]}'
        self.tsv_string = ''

    def test_getResponseContent(self):
        responseJSON = self.magma.getResponseContent(self.json_string)
        self.assertTrue(isinstance(responseJSON, dict))
        self.assertEqual(responseJSON['projects'][0]['project_name'], 'mvir1')

    def test_magmaCall(self):
        payload = {
            "project_name": 'ipi',
            "model_name": "sample",
            "record_names": ['IPIBLAD013'],
            "attribute_names": ['flojo_file']
        }
        magmaResponse, magmaResponseHeaders = self.magma.magmaCall(payload)
        self.assertTrue(isinstance(magmaResponse, dict))
        self.assertTrue(isinstance(magmaResponseHeaders, dict))
        self.assertEqual(magmaResponse, 'stuff') # TODO get actual data
        self.assertEqual(magmaResponseHeaders, 'headers_stuff') # TODO get actual data






    
