import requests
from requests import Session
from requests.auth import AuthBase
from requests.adapters import HTTPAdapter
from typing import Dict, Optional, List
from .auth import TokenAuth
from urllib3 import Retry
from urllib.parse import quote

def encode_path(*segments: str) -> str:
    return "/".join(quote(s) for s in segments)

class EtnaSession(Session):
    def __init__(self,
        auth: Optional[AuthBase],
        hooks: Dict = {}):
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

        for hook_name, hook_list in hooks.items():
            self.hooks[hook_name].extend(hook_list)

class EtnaClientBase:
    session: requests.Session
    hostname: str

    def __init__(self, session: EtnaSession, hostname: str):
        self.session = session
        self.hostname = hostname

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

