from typing import List
import pandas as pd
from functools import partial

from magby.Magma import *

_session = Session()

class Magby(object):
    def __init__(self, url: str, token: str):
        '''
        Client-facing class to interface with Magma (or Janus - for projects)
        :param url: str. Magma url
        :param token: str. Private token
        '''
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


    def retrieve(self,
                 projectName: str,
                 modelName: Union[List, str],
                 recordNames: Union[List, str]='all',
                 attributeNames: Union[List, str]='all',
                 dataType: str='meta',
                 filter='',
                 **kwargs) -> Union[pd.DataFrame, Dict]:
        '''
        Method to retrieve data from Magma
        :param projectName: str. Project name in Magma
        :param modelName: List[str] or str. Name(s) of model(s) in Magma
        :param recordNames: List[str] or str. Name(s) of record(s) in Magma
        :param attributeNames: List[str] or str. Name(s) of attribute(s) in Magma
        :param dataType: str. Type of data output. Allowed:
                        'meta' - metadata. Output - pandas.DataFrame
                        'json' - works for any type of data in Magma. Output - Dictionary
                        'mtx' - matrix of floats. Suitable for gene counts or  similar. Corresponds to 'matrix'
                                in Magma data types. Output - pandas.DataFrame
        :param filter: str. Filter to pass to Magma query. See magma docs
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: Data from magma in the form of pandas.DataFrame or Dict, depending on the dataType argument
        '''
        attributeNames = self._isAll(attributeNames)
        if (dataType == 'mtx') & (len(attributeNames) > 1):
            raise MagmaError('Magby.retrieve(): retrieval of matrix data is limited to a single attribute at a time')
        typeSelection = self._selectFormat(dataType)
        payload = self._constructPayload(projectName, modelName, recordNames, attributeNames,
                                         format=typeSelection[0], filter=filter)
        magma = self._recordMagmaObj(endpoint='retrieve', fmt=typeSelection[0], **kwargs)
        content, _ = self._call_api(payload, magma)
        if dataType == 'mtx':
            return self._wrapReturn(typeSelection[1], modelName=modelName, attributeName=attributeNames[0])(content)
        else:
            return self._wrapReturn(typeSelection[1])(content)


    def query(self,
              projectName: str,
              queryTerms: List,
              format: str = "json",
              **kwargs) -> Dict:
        '''
        Performs an expressive query to Magma.
        :param projectName: str. Project name in Magma
        :param queryTerms: List[str]. Query terms to magma. See Magma documentation
        :param format: str. Adjusts the format of the returned data if set to "tsv". See Magma documentation
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: Dict
        '''
        typeSelection = self._selectFormat('json')
        payload = {
            "project_name": projectName,
            "query": queryTerms,
        }
        if format=="tsv":
            typeSelection = self._selectFormat('meta')
            payload['format']="tsv"
        magma = self._recordMagmaObj(endpoint='query', fmt=typeSelection[0], **kwargs)
        content, _ = self._call_api(payload, magma)

        return self._wrapReturn(typeSelection[1])(content)


    def getProjects(self, **kwargs) -> List:
        '''
        Gets available projects in Magma for a given user (token)
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: List of Dict with abbreviated project names, full project names and descriptions
        '''
        janusUrl = self._janusUrl()
        magma = self._recordMagmaObj('retrieve', 'json', **kwargs)
        janusProjects = magma._janusProjectsCall(janusUrl)
        return janusProjects['projects']



    ### PRIVATE METHODS
    def _constructPayload(self,
                          projectName: str,
                          modelName: Union[List, str],
                          recordNames: Union[List, str],
                          attributeNames: Union[List, str],
                          format: str,
                          filter: str='') -> Dict:
        payload = {
            "project_name": projectName,
            "model_name": modelName,
            "record_names": self._isAll(recordNames),
            "attribute_names": self._isAll(attributeNames),
            "format": format,
            "filter": filter
        }
        return payload

    def _isAll(self, content: Union[List, str]) -> Union[List, str]:
        if content == 'all':
            return content
        else:
            return self._encapsulateToList(content)

    @staticmethod
    def _encapsulateToList(content: Union[List, str]) -> List:
        if isinstance(content, list):
            return content
        if isinstance(content, str):
            return [content]

    def _janusUrl(self) -> str:
        janusUrl = self._url.replace('magma', 'janus')
        return ('/'.join([janusUrl, 'projects']))


    def _recordMagmaObj(self, endpoint: str, fmt: str= 'json', session=_session) -> Magma:
        return Magma(url=self._url, token=self._token, endpoint=endpoint,
                     fmt=fmt, session=session)


    def _call_api(self, payload: Dict, magma: Magma):
        return magma.magmaCall(payload)

    @staticmethod
    def _selectFormat(dataType: str) -> List:
        '''
        Selector for types of outputs to retrieve from magma
        :param dataType: str
        :return: List where [0] is magma output format, [1] is a desired Magby output
        '''
        selection = {
            'meta': ['tsv', 'meta'],
            'json': ['json', 'json'],
            'mtx': ['json', 'mtx']
        }
        return selection[dataType]

    def _wrapReturn(self, fmt: str, **kwargs):
        switchReturns = {
            'meta': partial(pd.read_csv, sep='\t'),
            'json': dict,
            'mtx': partial(self._walkMatrix, **kwargs)
        }
        return switchReturns.get(fmt)

    @staticmethod
    def _walkMatrix(response: Dict, modelName: str, attributeName: str):
        observations = {x: response['models'][modelName]['documents'][x][attributeName] for x in response['models'][modelName]['documents']}
        features = response['models'][modelName]['template']['attributes'][attributeName]['options']
        return pd.DataFrame(observations, index=features)











    # TODO

    def _getModels(self, payload) -> List:
        pass

    def _getIDs(self, payload) -> List:
        pass

    def _getAttributes(self, payload) -> List:
        pass


    # TODO functions


    def _validUrl(self, url: str) -> bool:
        pass    # TODO implement URL validation



