from archimedes.functions.environment import token, magma_host, project_name
from archimedes.functions.magby import Magby
from archimedes.functions.list import flatten

def connect():
    return Magby.Magby(
        url=magma_host,
        token=token
    )

def question(magma, question, strip_identifiers=True):
    query_result = magma.query(project_name, queryTerms=question)

    if not 'answer' in query_result:
        raise Exception('No answer to magma query with elements: '+ ','.join(flatten(question)))

    return [ v[1] for v in query_result['answer'] ] if strip_identifiers else query_result['answer']
