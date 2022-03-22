from datetime import timedelta, datetime
from typing import Iterable, Tuple, Dict

from airflow import DAG
from airflow.decorators import task
from airflow.executors.debug_executor import DebugExecutor
from airflow.models import DagBag
from prometheus_client.samples import Sample

from etna.metrics.rollup_metrics import Rollup, RollupMetricsCollector
from etna import rollup, rollup_dag, system_dag


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
        ignore_first_depends_on_past=True,
    )

class TestRollup(Rollup):
    def __init__(self, value: int):
        self.value = value
        self.options = [
            (dict(a='aaaa', b='bbbb', c='cccc'), value)
        ]

    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        return self.options

def test_rollup_metrics():
    @rollup_dag(timedelta(minutes=1, seconds=60), ['a', 'b', 'c'])
    def test_rollup_dag():
        @task
        def do_thing():
            return TestRollup(1)

        @rollup
        def concat_test_rollups(a: TestRollup, b: TestRollup) -> TestRollup:
            return TestRollup(a.value + b.value)

        concat_test_rollups(do_thing())


    test_rollup_dag: DAG
    dag_bag = DagBag(read_dags_from_db=True)
    dag_bag.bag_dag(test_rollup_dag, test_rollup_dag)
    dag_bag.sync_to_db()

    run_dag(test_rollup_dag, datetime(2022, 1, 1, 1, 1), datetime(2022, 1, 1, 1, 10))

    samples = list(
        s
        for metric in RollupMetricsCollector().collect()
        for s in metric.samples
    )

    assert samples == [Sample('test_rollup_dag_rollup_test_rollup_dag_do_thing', labels=dict(a='aaaa', b='bbbb', c='cccc'), value=7, timestamp=None, exemplar=None)]