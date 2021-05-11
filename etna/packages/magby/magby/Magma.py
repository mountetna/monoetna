from io import StringIO
from typing import Dict, Tuple, Union
import json
from requests import Session

_session = Session()

class MagmaError(Exception):
    """
    errors corresponding to misuse of Magma
    """


class Magma(object):
    def __init__(self, url: str,
                 token: str,
                 endpoint: str,
                 fmt: str = 'json',
                 session: Session = _session) -> None:
        """
        Private class to interface with Magma
        :param url: Magma url
        :param token: Private Janus token
        :param endpoint: Magma (and janus) endpoints for HTTP requests
        :param fmt: Format of Magma output
        :param session: Instance of requests.Session. Needed to modify proxies and SSL certificate verification
        """

        allowedEndpoints = ['update', 'retrieve', 'query', 'update_model']
        if endpoint not in allowedEndpoints:
            raise MagmaError(f'Magma(): unknown endpoint {endpoint}. Must be one of the {allowedEndpoints}')
        self._url = '/'.join([url.strip('/'), endpoint])
        self._token = token
        self._fmt = fmt
        self._session = session
        self._headers = {"Content-Type": "application/json",
                         "Authorization": f"Etna {self._token}"}

    def getResponseContent(self, response) -> Union[Dict, StringIO]:
        switch = {
            'json': self._parseJSON,
            'tsv': self._parseText
        }
        return switch[self._fmt](response)

    def _parseJSON(self, response) -> Dict:
        return json.loads(response.text, strict=False)

    def _parseText(self, response) -> StringIO:
        return StringIO(response.text)

    def magmaCall(self, payload: Dict) -> Tuple:
        """
        Makes an API call to Magma and returns data and current headers
        :param payload: Dict. Payload
        :return: Tuple[Dict, Dict]. Dictionary of content, dictionary of response headers
        """
        response = self._session.post(
            self._url,
            data=json.dumps(payload),
            headers=self._headers)
        content = self.getResponseContent(response)
        return content, response.headers

    def _janusProjectsCall(self, url: str) -> Dict:
        """
        Get a projects object form Janus for a given user/token
        :param url: Janus/projects url
        :return: Dictionary of projects
        """
        response = self._session.get(
            url,
            headers=self._headers)
        content = self.getResponseContent(response)
        return content
