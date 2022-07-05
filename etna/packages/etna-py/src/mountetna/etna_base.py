import requests
from requests import Session, HTTPError
from requests.auth import AuthBase
from requests.adapters import HTTPAdapter
from typing import Dict, Optional, List, Type
from .auth import TokenAuth
from urllib3 import Retry
from urllib.parse import quote

def encode_path(*segments: str) -> str:
    return "/".join(quote(s) for s in segments)

class EtnaSession(Session):
    @staticmethod
    def assert_status(response: requests.Response, *args, **kwds):
        if 200 <= response.status_code < 300:
            return

        error_message = response.reason

        try:
            err_json = response.json()
            errors = err_json.get("errors", [])
            error = err_json.get("error")

            if error:
                errors.append(error)

            error_message = ", ".join(errors)
            error_message = (
                f"Request failed with {response.status_code}: {error_message}"
            )
        except (TypeError, ValueError, json.JSONDecodeError):
            pass

        raise HTTPError(error_message, response=response)

    def __init__(self, auth: Optional[AuthBase]):
        super().__init__()
        self.mount(
            "https://",
            HTTPAdapter(
                max_retries=Retry(
                    total=5,
                    backoff_factor=1,
                    connect=5,
                    read=3,
                    status_forcelist=[502, 503, 504],
                )
            )
        )
        if auth:
            self.auth = auth

        self.hooks['response'].extend([EtnaSession.assert_status])

class EtnaClientBase:
    auth: AuthBase
    session: requests.Session
    hostname: str

    def __init__(self, auth: Optional[AuthBase], hostname: str, session: Optional[Type[requests.Session]] = None):
        self.auth = auth
        self.session = (session or EtnaSession)(auth=auth)
        self.hostname = hostname.replace('https://', '')

    @property
    def auth_token(self) -> Optional[bytes]:
        if isinstance(self.session.auth, TokenAuth):
            return self.session.auth.token

    def prepare_url(self, *path: str):
        return f"https://{self.hostname}/{encode_path(*path)}"

    def prepare_url_unsafe(self, path: str):
        return f"https://{self.hostname}/{path}"

    def get_project_scope(self) -> Optional[str]:
        if isinstance(self.session.auth, TokenAuth):
            return self.session.auth.project_scope
        return None

