from unittest import TestCase
import configparser
from ..magby.Magby import *
from ..tests.testUtils import *

conf = configparser.ConfigParser()
conf.read('magby/tests/proxyConfig.ini')    # Config is in .gitignore



class TestMagby(TestCase):
    def setUp(self) -> None:
        self.proxy = json.loads(conf['DEFAULT'].get('proxy'))
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
                                               attributeNames=["patient"])
        self.assertTrue(isinstance(payload, dict))


    def test_getProjects(self):
        projects = self.magby.getProjects(proxies=self.proxy, verify=False)
        self.assertEqual(projects[2]['project_name'], 'ipi')
        self.assertTrue(isinstance(projects, list))
        self.assertTrue(len(projects) > 0)

    def test_retrieve(self):
        out = self.magby.retrieve(projectName="ipi", modelName='sample',
                                  recordNames='all', attributeNames=["patient"], dataType='df',
                                  proxies=self.proxy, verify=False)
        self.assertTrue(isinstance(out, pd.DataFrame))
        self.assertEqual(out.shape, (771,2))




