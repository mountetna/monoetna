def __getattr__(name):
    from .dags.decorators import system_dag
    from .operators import run_on_docker
    from .etls.decorators import metis_etl
    from .etls.metis import filter_by_record_directory, filter_by_exists_in_timur, list_contents_of_matches, link, MetisEtlHelpers
    from .xcom.etna_xcom import pickled
    from .hooks.hook_helpers import get_project_name, get_etna_hook, get_git_hook, get_project_slack_hook
    from .hooks.etna import UpdateRequest
    return locals()[name]


# trick the static checker.
if globals().get('notathing', False):
    from .dags.decorators import system_dag
    from .operators import run_on_docker
    from .etls.decorators import metis_etl
    from .etls.metis import filter_by_record_directory, filter_by_exists_in_timur, list_contents_of_matches, link, MetisEtlHelpers
    from .xcom.etna_xcom import pickled
    from .hooks.hook_helpers import get_project_name, get_etna_hook, get_git_hook, get_project_slack_hook
    from .hooks.etna import UpdateRequest
