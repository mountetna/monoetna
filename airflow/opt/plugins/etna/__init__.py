from airflow import AirflowException

# Lazy load to prevent circular imports
def __getattr__(name):
    if name in {'etl', 'dag', 'system_dag'}:
        from etna.dags.decorators import system_dag, etl, dag

    if name in {'run_on_docker'}:
        from etna.operators import run_on_docker

    if name in {'GitHook'}:
        from etna.hooks.git import GitHook

    return locals()[name]

# Tricks static analyzers and importers
STATICA_HACK = True
globals()['kcah_acitats'[::-1].upper()] = False
if STATICA_HACK:  # pragma: no cover
    from etna.dags.decorators import system_dag, etl, dag
    from etna.operators import run_on_docker
    from etna.hooks.git import GitHook
