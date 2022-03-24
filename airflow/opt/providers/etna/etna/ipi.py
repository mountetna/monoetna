from datetime import timedelta
from typing import List, Dict, Iterable, Tuple

from airflow.decorators import task
from airflow.operators.python import get_current_context
from airflow.utils.task_group import TaskGroup
from etna.etls.etl_task_batching import get_batch_range

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers, get_project_name, rollup_dag, rollup,
)
import re

from etna.hooks.etna import EtnaHook, Folder, File

from etna.metrics.rollup_metrics import Rollup

completion_metric_labels = ['attribute_name', 'attribute_type', 'model_name', 'project_name']

class MetisDataRollup(Rollup):
    project_name: str
    bucket_name: str
    folders: List[Folder]
    files: List[File]

    def __init__(self, project_name: str, bucket_name: str):
        self.project_name = project_name
        self.bucket_name = bucket_name
        self.folders = []
        self.files = []

    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        for file in files:
            pass

    def concat(self, other: "MetisDataRollup") -> "MetisDataRollup":
        if self.project_name != other.project_name or self.bucket_name != other.bucket_name:
            raise ValueError('Cannot rollup metis data from different buckets')

        result = MetisDataRollup(self.project_name, self.bucket_name)

        return result

class MagmaCompletionRollup(Rollup):
    project_name: str
    all_record_ids: Dict[str, set]
    attribute_types: Dict[str, Dict[str, str]]
    complete_attributes: Dict[str, Dict[str, set]]
    sampled_record_ids: Dict[str, set]

    def __init__(self, project_name: str):
        self.project_name = project_name
        self.all_record_ids = {}
        self.attribute_types = {}
        self.complete_attributes = {}
        self.sampled_record_ids = {}

    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        for model_name, attributes in self.attribute_types.items():
            for attribute_name in attributes.keys():
                row_ids = self.complete_attribute_record_ids(model_name, attribute_name) & self.all_record_ids.setdefault(model_name, set())
                yield dict(
                    model_name=model_name,
                    attribute_name=attribute_name,
                    attribute_type=self.attribute_type(model_name, attribute_name, 'unknown'),
                ), len(row_ids)

    def sample_record_id(self, model_name: str, record_name: str):
        self.sampled_record_ids.setdefault(model_name, set()).add(record_name)

    def complete_attribute_record_ids(self, model_name: str, attribute_name: str) -> set:
        return self.complete_attributes.setdefault(model_name, {}).setdefault(attribute_name, set())

    def incomplete_attribute_record_ids(self, model_name: str, attribute_name: str) -> set:
        return self.sampled_record_ids.get(model_name, set()) - self.complete_attribute_record_ids(model_name, attribute_name)

    def attribute_type(self, model_name: str, attribute_name: str, default_type: str) -> str:
        return self.attribute_types.setdefault(model_name, {}).setdefault(attribute_name, default_type)

    def concat(self, other: "MagmaCompletionRollup") -> "MagmaCompletionRollup":
        result = MagmaCompletionRollup(self.project_name)
        if self.project_name != other.project_name:
            raise ValueError(f"Projects {self.project_name} and {other.project_name} cannot be concatenated")

        result.all_record_ids = other.all_record_ids
        result.attribute_types = other.attribute_types
        result.complete_attributes = self.complete_attributes.copy()

        for model_name, model_attribute_types in self.attribute_types.items():
            for attribute_name, attribute_type in model_attribute_types.items():
                for id in other.complete_attribute_record_ids(model_name, attribute_name):
                    result.complete_attribute_record_ids(model_name, attribute_name).add(id)

                for id in other.incomplete_attribute_record_ids(model_name, attribute_name):
                    result.complete_attribute_record_ids(model_name, attribute_name).remove(id)

        return result

def collect_metis_metrics(project_name: str, bucket_name: str):
    etna_hook = EtnaHook.for_project(project_name)

    def process_project_metis_data():
        @task
        def gather_updates():
            rollup = MetisDataRollup(project_name, bucket_name)

            with etna_hook.metis(project_name) as metis:
                start, end = get_batch_range(get_current_context())
                files, _ = metis.tail(project_name, bucket_name, 'files', batch_start=start, batch_end=end)
                _, folders = metis.tail(project_name, bucket_name, 'folders', batch_start=start, batch_end=end)
                rollup.files.extend(files)
                rollup.folders.extend(folders)

            return rollup

        @rollup
        def concat_metis_rollups(a: MetisDataRollup, b: MetisDataRollup) -> MetisDataRollup:
            return a.concat(b)

    process_project_metis_data.__name__ = f"process_{project_name}_metis_data"
    rollup_dag(timedelta(minutes=5), completion_metric_labels)(process_project_metis_data)

def collect_completion_metrics(project_name: str):
    etna_hook = EtnaHook.for_project(project_name)

    def process_project_magma_completion():
        @task
        def gather_updates():
            rollup = MagmaCompletionRollup(project_name)

            with etna_hook.magma(project_name=project_name) as magma:
                models = magma.retrieve(
                    project_name=project_name,
                    model_name="all",
                    hide_templates=False,
                ).models

                for model_name, model in models.items():
                    for attribute_name, attribute in model.template.attributes.items():
                        rollup.attribute_type(model_name, attribute_name, attribute.attribute_type)

                    document_model = magma.retrieve(
                        project_name=project_name,
                        model_name=model_name,
                        attribute_names="identifier"
                    ).models.get(model_name)

                    if not document_model:
                        continue

                    rollup.all_record_ids.setdefault(model_name, set(document_model.documents.keys()))

                    start, end = get_batch_range(get_current_context())
                    model_updates = magma.retrieve(
                        project_name=project_name,
                        model_name=model_name,
                        attribute_names='all',
                        record_names='all',
                        order='updated_at',
                        filter=[
                            f"updated_at>={start.isoformat(sep='seconds')}"
                            f"updated_at<={start.isoformat(sep='seconds')}"
                        ]
                    ).models.get(model_name)

                    if not model_updates:
                        continue

                    for record_name, record in model_updates.documents.items():
                        rollup.sample_record_id(model_name, record_name)
                        for attribute_name in model.template.attributes.keys():
                            completed = rollup.complete_attribute_record_ids(model_name, attribute_name)
                            if record.get(attribute_name) is not None:
                                completed.add(record_name)

        @rollup
        def concat_magma_completion_rollups(a: MagmaCompletionRollup, b: MagmaCompletionRollup) -> MagmaCompletionRollup:
            return a.concat(b)

    process_project_magma_completion.__name__ = f"process_{project_name}_completion_metrics"
    rollup_dag(timedelta(minutes=5), completion_metric_labels)(process_project_magma_completion)