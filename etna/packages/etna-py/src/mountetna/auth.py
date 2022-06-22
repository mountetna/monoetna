import requests
from typing import Optional
from requests.auth import AuthBase
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.asymmetric import padding
from requests.auth import AuthBase
import base64

class TokenAuth(AuthBase):
    token: bytes
    project_scope: Optional[str]

    def __init__(self, token: bytes, project_scope: Optional[str] = None):
        self.token = token
        self.project_scope = project_scope

    def __call__(self, r: requests.Request) -> requests.Request:
        r.headers["Authorization"] = f'Etna {self.token.decode("ascii")}'
        return r

class SigAuth(AuthBase):
    nonce: bytes
    email: str
    private_key: bytes

    def __init__(self, private_key_file: str, email: str, nonce: bytes):
        self.private_key = open(private_key_file, "rb").read()
        self.nonce = nonce
        self.email = email

    def __call__(self, r: requests.Request) -> requests.Request:
        key = serialization.load_pem_private_key(self.private_key, None)
        txt_to_sign: bytes = b".".join(
            [self.nonce, base64.b64encode(self.email.encode("ascii"))]
        )
        sig = base64.b64encode(
            key.sign(
                txt_to_sign,
                padding.PKCS1v15(
                    # If in the future the server changes to PSS
                    # mgf=padding.MGF1(hashes.SHA256()),
                    # salt_length=padding.PSS.MAX_LENGTH
                ),
                hashes.SHA256(),
            )
        )
        together: bytes = b".".join([txt_to_sign, sig])
        r.headers["Authorization"] = f'Signed-Nonce {together.decode("ascii")}'
        return r

