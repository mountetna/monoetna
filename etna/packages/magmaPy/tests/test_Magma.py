from unittest import TestCase
from ..magmaPy.Magma import *


class TestMagma(TestCase):
    def setUp(self) -> None:
        self.magma = Magma('some_url', 'my_token')
        self.json_string = '{"projects":[{"project_name":"mvir1","project_name_full":"COMET","project_description":null},' \
                           '{"project_name":"dscolab","project_name_full":"Data Science CoLab","project_description":null}]}'

    def test_url(self):
        self.assertEqual(self.magma.url, 'some_url')
        self.magma.url = 'changed_url'
        self.assertNotEqual(self.magma.url, 'some_url')
        self.assertEqual(self.magma.url, 'changed_url')
        
    def test_token(self):
        self.assertEqual(self.magma.token, 'my_token')
        self.magma.token = 'changed_token'
        self.assertNotEqual(self.magma.token, 'my_token')
        self.assertEqual(self.magma.token, 'changed_token')

    def test_extractContent(self):
        responseContent = self.magma.getContent(self.json_string)
        self.assertEqual(responseContent['projects'][0]['project_name'], 'mvir1')

    def test_magmaCall(self):
        pass




    
