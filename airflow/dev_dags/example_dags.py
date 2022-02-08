from datetime import timedelta
from etna import system_dag, run_on_docker

@system_dag(timedelta(days=1))
def system_dag_example():
    run_on_docker('test_container_operator', 'polyphemus_app',
                     ["/app/bin/polyphemus", "etl", "ipi_watch_files_etl", "run"])