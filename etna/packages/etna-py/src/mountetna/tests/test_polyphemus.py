from .. import Polyphemus, EtnaSession, TokenAuth
from ..polyphemus import ListEtlConfigsResponse
import pytest
import responses
from typing import Dict
import json
from datetime import datetime, timedelta

@responses.activate
def test_polyphemus_etl_config_should_run():
    ''' it should run in the indicated dates'''
    ran_at = datetime.now() - timedelta(1)
    start = datetime.now()
    end = datetime.now() + timedelta(1)

    configs = {
        'configs': [
            {
                'project_name': 'labors',
                'ran_at': (ran_at).isoformat(),
                'run_interval': 600
            }
        ]
    }
    responses.add(
        responses.POST,
        "https://polyphemus.test/api/etl/configs",
        body=json.dumps(configs),
        status=200
    )
    client = Polyphemus(
        auth=TokenAuth(token='token'),
        hostname='polyphemus.test'
    )
    response = client.list_all_etl_configs()
    config = response.configs[0]

    assert config.should_run(start, end)
    assert not config.should_run(ran_at, end)
