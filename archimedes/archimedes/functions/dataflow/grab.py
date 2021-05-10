from requests import Session
from archimedes.functions.environment import token

def curl_data(url, token = token, verify = True):
    get = Session().get(
        url,
        headers={
            "Content-Type" : "application/json",
            "Authorization" : f"Etna {token}"},
        verify=verify)
    get.raise_for_status()
    return get