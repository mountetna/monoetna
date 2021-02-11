from unittest import TestCase
import configparser
from ..magby.Magby import *
from ..tests.testUtils import *

conf = configparser.ConfigParser()
conf.read('magby/tests/proxyConfig.ini')    # Config is in .gitignore



class TestMagby(TestCase):
    def setUp(self) -> None:
        proxy = json.loads(conf['DEFAULT'].get('proxy'))
        self.magby = Magby(conf['DEFAULT'].get('url'),
                               conf['DEFAULT'].get('token'))


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

    def test_getProjects(self):
        projects = self.magby.getProjects()
        self.assertTrue('ipi' in projects)
        self.assertTrue(isinstance(projects, list))
        self.assertTrue(len(projects) > 0)

    @vcr.use_cassette()
    def test_retrieve(self):
        payload = self.magby.constructPayload()
        out = self.magby.retrieve(payload, 'df')
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (10,10))





