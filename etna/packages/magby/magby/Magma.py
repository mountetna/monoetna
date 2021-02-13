from io import StringIO
from typing import Dict, Tuple, Union
import json
from requests import RequestException, Session


MagmaError = RequestException

_session = Session()


class MagmaError(Exception):
    '''
    errors corresponding to misuse of Magma
    '''


class Magma(object):
    def __init__(self, url: str,
                 token: str,
                 endpoint: str,
                 fmt: str='json',
                 session: Session()=_session) -> None:
        allowedEndpoints = ['update', 'retrieve', 'query', 'update_model']
        if endpoint not in allowedEndpoints:
            raise MagmaError(f'Magma(): unknown endpoint {endpoint}. Must be one of the {allowedEndpoints}')
        self._url = '/'.join([url.strip('/'), endpoint])
        self._token = token
        self._fmt = fmt
        self._session = session
        self._headers = {"Content-Type" : "application/json",
                         "Authorization" : f"Etna {self._token}"}


    def getResponseContent(self, response) -> Union[Dict, StringIO]:
        '''
        An abstraction to
        :param response: magma API response
        :return: Dict
        '''
        switch = {
            'json': self._parseJSON,
            'tsv': self._parseText
        }

        return switch[self._fmt](response)


    def _parseJSON(self, response) -> Dict:
        return json.loads(response.text, strict=False)


    def _parseText(self, response) -> StringIO:
        return StringIO(response.text)


    def magmaCall(self, payload: Dict, **kwargs) -> Tuple:
        '''
        Make API call and return data and current headers
        :param payload: Dict. Payload
        :param **kwargs: passed to requests.Session.post()
        :return: Tuple[Dict, Dict]. Dictionary of content, dictionary of response headers
        '''
        response = self._session.post(self._url, data=json.dumps(payload), headers=self._headers, **kwargs)
        content = self.getResponseContent(response)
        return content, response.headers

    def _janusProjectsCall(self, url: str, **kwargs) -> Dict:
        '''
        Get a projects object form Janus for a given user/token
        :param url: Janus/projects url
        :param kwargs: kwargs passed to requests.Session() methods
        :return: Dictionary of projects
        '''
        response = self._session.get(url, headers=self._headers, **kwargs)
        content = self.getResponseContent(response)
        return content





