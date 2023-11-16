from .. import Magma, TokenAuth, Model, UpdateRequest
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
def test_magma_update(project_template: Dict, labor_template: Dict):
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
				},
				"project": {
                    "template": project_template
                }
			}
        },
        status=200
    )
    client = Magma(
        auth=TokenAuth(token='token'),
        hostname='magma.test'
    )
    response = client.update(UpdateRequest(**{ "project_name":"ipi", "revisions": { "labor": { "The Nemean Lion" : { "country" : "Nemea" } } } }) )
    assert response.models['labor'].documents['The Nemean Lion']['country'] == 'Nemea'

@responses.activate
def test_magma_update_pages(labor_template: Dict):
    ''' it updates records in magma by pages, with dry_run and autolink passed through'''
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
    response = client.update(UpdateRequest(**{ "project_name":"ipi", "revisions": { "labor": {
        "The Nemean Lion" : { "country" : "Nemea" },
        "The Lernean Hydra" : { "country" : "Lerna" },
        "The Ceryneian Hind" : { "country" : "Ceryneia" },
        "The Erymanthian Boar" : { "country" : "Erymanthia" }
    } },
    "dry_run": True, "autolink": True}), page_size=2 )

    assert responses.assert_call_count("https://magma.test/update", 2) is True
    len_calls=len(responses.calls)
    last2=[len_calls-2, len_calls-1]
    assert all(['"dry_run":true' in str(responses.calls[i].request.body) for i in last2]) is True
    assert all(['"autolink":true' in str(responses.calls[i].request.body) for i in last2]) is True

    response = client.update(UpdateRequest(**{ "project_name":"ipi", "revisions": { "labor": {
        "The Nemean Lion" : { "country" : "Nemea" },
        "The Lernean Hydra" : { "country" : "Lerna" },
        "The Ceryneian Hind" : { "country" : "Ceryneia" },
        "The Erymanthian Boar" : { "country" : "Erymanthia" }
    } },
    "dry_run": False, "autolink": False}), page_size=2 )

    assert responses.assert_call_count("https://magma.test/update", 4) is True
    len_calls=len(responses.calls)
    last2=[len_calls-2, len_calls-1]
    assert all(['"dry_run":false' in str(responses.calls[i].request.body) for i in last2]) is True
    assert all(['"autolink":false' in str(responses.calls[i].request.body) for i in last2]) is True

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


@responses.activate
def test_magma_query_tsv():
    ''' it queries magma '''
    body = "identifier\tcolumn1\tcolumn2\nrecord1\t1\tblahblah\nrecord2\t3.4\tanything you want"
    responses.add(
        responses.POST,
        "https://magma.test/query",
        body=body,
        status=200
    )
    client = Magma(
        auth=TokenAuth(token='token'),
        hostname='magma.test'
    )

    response = client.query({ "project_name":"ipi", "query": [ 'labor', [ 'country', '::equals', 'Nemea' ], '::all', 'country' ], "format": "tsv"})
    assert response.to_json() == '{"identifier":{"0":"record1","1":"record2"},"column1":{"0":1.0,"1":3.4},"column2":{"0":"blahblah","1":"anything you want"}}'
