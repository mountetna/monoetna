from datetime import timedelta
from typing import List, Dict, Iterable, Tuple

from airflow.decorators import task
from airflow.operators.python import get_current_context
from airflow.utils.task_group import TaskGroup
from etna.etls.etl_task_batching import get_batch_range

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers, get_project_name, rollup_dag,
)
import re

from etna.etls.metis import MatchedRecordFolder
from etna.hooks.etna import Attribute, UpdateRequest, EtnaHook
import csv
import json
import os.path

from etna.metrics.rollup_metrics import Rollup

completion_metric_labels = ['attribute_name', 'attribute_type', 'model_name', 'project_name']

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

    def attribute_type(self, model_name: str, attribute_name: str, default_type: str) -> str:
        return self.attribute_types.setdefault(model_name, {}).setdefault(attribute_name, default_type)

    def concat(self, other: "MagmaCompletionRollup") -> "MagmaCompletionRollup":
        result = MagmaCompletionRollup(self.project_name)
        if self.project_name != other.project_name:
            raise ValueError(f"Projects {self.project_name} and {other.project_name} cannot be concatenated")

        return result


def measure_completion_metrics(project_name: str):
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

    process_project_magma_completion.__name__ = f"process_{project_name}_completion_metrics"
    rollup_dag(timedelta(minutes=5), completion_metric_labels)(process_project_magma_completion)