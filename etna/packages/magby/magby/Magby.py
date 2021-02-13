from typing import List, Dict
import json
import pandas as pd
import numpy as np
import os
import pickle
import pprint
from functools import partial

from ..magby.Magma import *

class Magby(object):
    def __init__(self, url, token):
        self._url = url.strip('/')
        self._token = token
        self._magma = None

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
        magma = self._formMagmaObj('retrieve')
        janusProjects = magma._janusProjectsCall(janusUrl, **kwargs)
        return janusProjects['projects']


    def _formMagmaObj(self, endpoint: str, fmt: str='json') -> Magma:
        return Magma(self._url, self._token, endpoint, fmt)


    def _call_api(self, payload: Dict, endpoint: str, fmt: str='json', **kwargs):
        magma = self._formMagmaObj(endpoint, fmt)
        return magma.magmaCall(payload, **kwargs)





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
        content, _ = self._call_api(payload, 'retrieve', fmt=typeSelection[0], **kwargs)
        switchReturns = {
            'meta': partial(pd.DataFrame, sep='\t'),
            'json': str,
            'mtx': np.array
        }
        return switchReturns[typeSelection[1]]


    def _selectFormat(self, dataType: str) -> List:
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



