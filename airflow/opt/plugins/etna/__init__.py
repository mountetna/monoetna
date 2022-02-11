from airflow import AirflowException

# Lazy load to prevent circular imports
def __getattr__(name):
    if name in {'etl', 'dag', 'system_dag'}:
        from etna.dags.decorators import system_dag, etl, dag

    if name in {'run_on_docker'}:
        from etna.operators import run_on_docker

    if name in {'GitHook'}:
        from etna.hooks.git import GitHook

    if name in {'pickled'}:
        from etna.xcom.etna_xcom import pickled

    if name in {'EtnaHook', 'Folder', 'File'}:
        from etna.hooks.etna import EtnaHook, Folder, File

    if name in {'filter_by_record_directory', 'MatchedAtRoot', 'filter_by_exists_in_timur', 'link', 'list_contents_of_matches'}:
        from etna.etls.metis import filter_by_record_directory, MatchedAtRoot, filter_by_exists_in_timur, link, \
            list_contents_of_matches

    if name in {'metis_etl'}:
        from etna.dags.decorators import metis_etl

    return locals()[name]

# Tricks static analyzers and importers
STATICA_HACK = True
globals()['kcah_acitats'[::-1].upper()] = False
if STATICA_HACK:  # pragma: no cover
    from etna.dags.decorators import system_dag, etl, dag
    from etna.operators import run_on_docker
    from etna.hooks.git import GitHook
    from etna.xcom.etna_xcom import pickled
    from etna.hooks.etna import EtnaHook, Folder, File
    from etna.etls.metis import filter_by_record_directory, MatchedAtRoot, filter_by_exists_in_timur, link, \
        list_contents_of_matches
    from etna.dags.decorators import metis_etl
