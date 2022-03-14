from datetime import timedelta, datetime
from typing import Iterable, Tuple, Dict

from airflow import DAG
from airflow.decorators import task
from airflow.executors.debug_executor import DebugExecutor
from prometheus_client.samples import Sample

from etna.dags.decorators import rollup_dag
from etna.metrics.rollup_metrics import Rollup, RollupMetricsCollector


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
        ignore_first_depends_on_past=True,
    )

class TestRollup(Rollup):
    def __init__(self):
        self.options = [
            (dict(a='aaaa', b='bbbb', c='cccc'), 159)
        ]

    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        return self.options

def test_rollup_metrics():
    @rollup_dag(timedelta(minutes=1), ['a', 'b', 'c'])
    def test_rollup_dag():
        @task
        def do_thing():
            return TestRollup()

        do_thing()

    test_rollup_dag: DAG

    run_dag(test_rollup_dag, datetime(2022, 1, 1, 1, 1), datetime(2022, 1, 1, 1, 2))

    for metric in RollupMetricsCollector().collect():
        assert metric.samples == [Sample('', dict(), 159, None)]