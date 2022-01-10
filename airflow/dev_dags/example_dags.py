from datetime import timedelta
from etna import run_in_container, system_dag


@system_dag(timedelta(days=1))
def system_dag_example():
    run_in_container('test_container_operator', 'polyphemus_app',
                     ["/app/bin/polyphemus", "etl", "ipi_watch_files_etl", "run"])