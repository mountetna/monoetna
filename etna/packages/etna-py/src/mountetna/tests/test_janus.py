from .. import Janus, EtnaSession, TokenAuth
import pytest
import responses
from typing import Dict

@responses.activate
def test_janus_generate_token():
    ''' it generates a task token via Janus '''
    responses.add(
        responses.POST,
        "https://janus.test/api/tokens/generate",
        body='task_token',
        status=200
    )
    client = Janus(
        auth=TokenAuth(token='token'),
        hostname='janus.test'
    )
    assert client.generate_token('labors') == 'task_token'
