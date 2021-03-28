from archimedes.functions.environment import token, magma_host, project_name, app_env
from archimedes.functions.magby import Magby

def connect():
    return Magby.Magby(
        url=magma_host,
        token=token,
        verify=(app_env == 'production')
    )

def question(magma, question, strip_identifiers=True):
    query_result = magma.query(project_name, queryTerms=question)

    if not 'answer' in query_result:
        raise 'No answer to magma query'

    return [ v[1] for v in query_result['answer'] ] if strip_identifiers else query_result['answer']
