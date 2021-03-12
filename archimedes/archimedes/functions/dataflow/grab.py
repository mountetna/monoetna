from requests import Session
from archimedes.functions.environment import token

def curl_data(url, token = token, verify = True):
    return Session().get(
        url,
        headers={
            "Content-Type" : "application/json",
            "Authorization" : f"Etna {token}"},
        verify=verify)