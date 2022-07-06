from .. import Metis, EtnaSession, TokenAuth
import pytest
import responses
from responses import GET, POST, add
from typing import Dict
from io import BytesIO

@responses.activate
def test_metis_upload_file():
    ''' it sends a file to metis '''
    add(GET, "https://metis.test/athena/list/plans/helmet", json={}, status=422)
    add(POST, "https://metis.test/athena/folder/create/plans/helmet", json={}, status=200)
    add(POST, "https://metis.test/authorize/upload", json={'url': 'upload/plans/helmet/helmet.txt'}, status=200)
    add(POST, "https://metis.test/upload/plans/helmet/helmet.txt", json={}, status=200)
    client = Metis(
        auth=TokenAuth(token='token'),
        hostname='metis.test'
    )
    helmet_file = BytesIO(b'''
      xXx
     xO|Ox
    ''')
    for upload in client.upload_file('athena', 'plans', 'helmet/helmet.txt', helmet_file):
        print('Uploading')
