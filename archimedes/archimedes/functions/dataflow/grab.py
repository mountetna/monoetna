from requests import Session
from archimedes.functions.environment import token

def curl_data(url, token = token):
    return Session().get(
        url,
        headers={
            "Content-Type" : "application/json",
            "Authorization" : f"Etna {token}"})