from .. import Magma, TokenAuth, Model
import pytest
import responses
from typing import Dict

@responses.activate
def test_magma_retrieve(labor_template: Dict):
    ''' it retrieves records from magma '''
    responses.add(
        responses.POST,
        "https://magma.test/retrieve",
        json={
			"models": {
				"labor": {
					"documents": {},
                    "template": labor_template
				}
			}
        },
        status=200
    )
    client = Magma(
        auth=TokenAuth(token='token'),
        hostname='magma.test'
    )
    response = client.retrieve('ipi')

    assert 'labor' in response.models
    assert type(response.models['labor']) == Model

@responses.activate
def test_magma_update(labor_template: Dict):
    ''' it updates records in magma '''
    responses.add(
        responses.POST,
        "https://magma.test/update",
        json={
			"models": {
				"labor": {
					"documents": {
                        "The Nemean Lion" : {
                            "country" : "Nemea"
                        }
                    },
                    "template": labor_template
				}
			}
        },
        status=200
    )
    client = Magma(
        auth=TokenAuth(token='token'),
        hostname='magma.test'
    )
    response = client.update({ "project_name":"ipi", "revisions": { "labor": { "The Nemean Lion" : { "country" : "Nemea" } } } } )
    assert response.models['labor'].documents['The Nemean Lion']['country'] == 'Nemea'

@responses.activate
def test_magma_query():
    ''' it queries magma '''
    responses.add(
        responses.POST,
        "https://magma.test/query",
        json={
            'answer': [ [ 'The Nemean Lion', 'Nemea' ] ],
            'format': [ 'labor#country' ]
        },
        status=200
    )
    client = Magma(
        auth=TokenAuth(token='token'),
        hostname='magma.test'
    )
    response = client.query({ "project_name":"ipi", "query": [ 'labor', [ 'country', '::equals', 'Nemea' ], '::all', 'country' ]})
    assert response.answer == [ [ 'The Nemean Lion', 'Nemea' ] ]
