from typing import List, Dict
import json
import pandas as pd
import os
import pickle
import pprint

from ..magby.Magma import *

class Magby(object):
    def __init__(self, url, token):
        self._url = url.strip('/')
        self._token = token

    @property
    def url(self) -> str:
        return self._url

    @url.setter
    def url(self, newUrl: str) -> None:
        self._url = newUrl

    @property
    def token(self) -> str:
        return self._token

    @token.setter
    def token(self, newToken: str) -> None:
        self._token = newToken

    @token.deleter
    def token(self) -> None:
        del self._token


    def _call_api(self, payload: Dict, fmt: str):
        magma = Magma(self.url, self.token, fmt)
        return magma.magmaCall(payload)



    # TODO

    def retrieve(self, payload: Dict, dataType: str, query: List=None) -> Union[pd.DataFrame, Tuple]:
        pass

    def constructPayload(self) -> Dict:
        pass

    def _getModels(self, payload) -> List:
        pass

    def _getIDs(self, payload) -> List:
        pass

    def _getAttributes(self, payload) -> List:
        pass








    # TODO functions
    def getProjects(self) -> List[str]:
        pass    # TODO get project names
    def _janusProjects(self) -> str:
        janusUrl = self._url.replace('magma', 'janus')
        return ('/'.join([janusUrl, 'projects']))
    def _validUrl(self, url: str) -> bool:
        pass    # TODO implement URL validation



