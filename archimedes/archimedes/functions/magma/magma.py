from archimedes.functions.environment import token, magma_host, project_name
from archimedes.functions.etna import Magma
from archimedes.functions.list import flatten

def connect():
    return Magma(
        session=EtnaSession(auth=TokenAuth(token=bytes(token))),
        hostname=magma_host
    )

def question(magma, question, strip_identifiers=True):
    query_result = magma.query({
        'project_name': project_name,
        'query': question
        })

    if not query_result.answer:
        raise Exception('No answer to magma query with elements: '+ ','.join(flatten(question)))

    return [ v[1] for v in query_result.answer ] if strip_identifiers else query_result.answer

def query_tsv(magma, project_name, queryTerms, user_columns=[], expand_matrices=False, transpose=False):
    query_result = magma.query({
        'project_name': project_name,
        'query': queryTerms,
        'format': 'tsv',
        'user_columns': user_columns,
        'expand_matrices': expand_matrices,
        'transpose': transpose)
    return query_result
