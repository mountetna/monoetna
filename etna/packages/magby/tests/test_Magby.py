from unittest import TestCase
from ..magby.Magby import *

class TestMagby(TestCase):
    def setUp(self) -> None:
        self.magby = Magby('some_url', 'my_token')
        self.json_string = '{"projects":[{"project_name":"mvir1","project_name_full":"COMET","project_description":null},' \
                           '{"project_name":"dscolab","project_name_full":"Data Science CoLab","project_description":null}]}'

    def test_url(self):
        self.assertEqual(self.magby.url, 'some_url')
        self.magby.url = 'changed_url'
        self.assertNotEqual(self.magby.url, 'some_url')
        self.assertEqual(self.magby.url, 'changed_url')

    def test_token(self):
        self.assertEqual(self.magby.token, 'my_token')
        self.magby.token = 'changed_token'
        self.assertNotEqual(self.magby.token, 'my_token')
        self.assertEqual(self.magby.token, 'changed_token')


