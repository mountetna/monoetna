import re
from datetime import datetime
from typing import Generator, Iterable, Tuple, Dict

from airflow.models import DagBag, XCom, BaseOperator
from prometheus_client import REGISTRY, Metric
from prometheus_client.metrics_core import GaugeMetricFamily
from prometheus_client.samples import Sample
from pytz import utc


rollup_label_metric_re = re.compile(r'label:(.*)')

# magma_record_attributes[model_name, restricted, attribute_name, attribute_type] = count

class Rollup:
    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        return []

class RollupMetricsCollector:
    def describe(self):
        return []

    def collect(self) -> Generator[Metric, None, None]:
        bag = DagBag(read_dags_from_db=True)
        bag.collect_dags_from_db()
        for dag_name, dag in bag.dags.items():
            if 'rollup' not in dag.tags:
                continue

            for task_id in dag.task_ids:
                task = dag.get_task(task_id)
                if task.task_type != 'RollupXcomOperator':
                    continue

                metric = GaugeMetricFamily(
                    f"{dag_name}_{task_id}",
                    dag.description or f"Quantity pulled from {dag_name}.{task_id}",
                    labels=[
                        match.group(1) for match in
                        (rollup_label_metric_re.match(t) for t in dag.tags)
                        if match is not None
                    ],
                )

                now = datetime.now().replace(tzinfo=utc)
                rollup_value = XCom.get_one(dag_id=dag_name, task_id=task_id, include_prior_dates=True, execution_date=now)

                if rollup_value is not None:
                    for labels, value in rollup_value.measure():
                        metric.samples.append(Sample(
                            metric.name, labels, value, None
                        ))

                yield metric

def register_rollup_metrics():
    REGISTRY.register(RollupMetricsCollector())