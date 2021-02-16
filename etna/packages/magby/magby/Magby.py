from typing import List, Dict
import json
import pandas as pd
import numpy as np
import os
import pickle
import pprint
from functools import partial

from ..magby.Magma import *

_session = Session()


class Magby(object):
    def __init__(self, url: str, token: str, session=_session):
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

    @staticmethod
    def _constructPayload(projectName: str,
                          modelName: Union[List, str],
                          recordNames: Union[List, str],
                          attributeNames: Union[List, str],
                          format: str,
                          filter: str=None) -> Dict:
        payload = {
            "project_name": projectName,
            "model_name": modelName,
            "record_names": recordNames,
            "attribute_names": attributeNames,
            "format": format,
            "filter": filter
        }
        return payload


    def _janusUrl(self) -> str:
        janusUrl = self._url.replace('magma', 'janus')
        return ('/'.join([janusUrl, 'projects']))


    def getProjects(self, **kwargs) -> List:
        janusUrl = self._janusUrl()
        magma = self._recordMagmaObj('retrieve')
        janusProjects = magma._janusProjectsCall(janusUrl, **kwargs)
        return janusProjects['projects']


    def _recordMagmaObj(self, endpoint: str, fmt: str= 'json', session=_session) -> Magma:
        return Magma(url=self._url, token=self._token, endpoint=endpoint,
                     fmt=fmt, session=session)


    def _call_api(self, payload: Dict, magma: Magma, **kwargs):
        return magma.magmaCall(payload, **kwargs)

    @staticmethod
    def _selectFormat(dataType: str) -> List:
        '''
        Selector for types of outputs to retrieve from magma
        :param dataType: str
        :return: List where [0] is magma output format, [1] is a desired Magby output
        '''
        selection = {
            'meta': ['tsv', 'df'],
            'json': ['json', 'json'],
            'mtx': ['json', 'mtx']
        }
        return selection[dataType]


###### NOT TESTED




    def retrieve(self,
                 projectName: str,
                 modelName: Union[List, str],
                 recordNames: Union[List, str]='all',
                 attributeNames: Union[List, str]='all',
                 dataType: str='df',
                 **kwargs) -> Union[pd.DataFrame, List, Dict]:

        typeSelection = self._selectFormat(dataType)
        payload = self._constructPayload(projectName, modelName, recordNames, attributeNames, format=typeSelection[0])
        magma = self._recordMagmaObj(endpoint='retrieve', fmt=typeSelection[0], **kwargs)
        content, _ = self._call_api(payload, magma, **kwargs)
        switchReturns = {
            'meta': partial(pd.DataFrame, sep='\t'),
            'json': str,
            'mtx': np.array
        }
        return switchReturns[typeSelection[1]]








    # TODO




    def _getModels(self, payload) -> List:
        pass

    def _getIDs(self, payload) -> List:
        pass

    def _getAttributes(self, payload) -> List:
        pass

    def query(self, queryTerms: List):
        pass

    # TODO functions


    def _validUrl(self, url: str) -> bool:
        pass    # TODO implement URL validation



