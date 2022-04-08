from datetime import datetime
from inspect import isgenerator
from typing import Dict, Optional, List
import re
from urllib.parse import quote
import logging

from airflow import DAG
from airflow.operators.python import get_current_context
from requests.adapters import HTTPAdapter
from serde import serialize, deserialize

from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from airflow.models.xcom_arg import XComArg

from etna.etls.etl_task_batching import get_batch_range

import cached_property
import dateutil
import requests
from airflow.exceptions import AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.asymmetric import padding
from requests import HTTPError, Session
from requests.auth import AuthBase
from serde.json import from_json, to_json
from urllib3 import Retry

from etna.dags.project_name import get_project_name
from etna.hooks.keys import prepared_key_from
from etna.utils.iterables import batch_iterable
from etna.utils.multipart import encode_as_multipart
from etna.utils.streaming import iterable_to_stream
from etna.hooks.box import BoxHook, BoxFile, Box


class BoxEtlHelpers:
    hook: BoxHook

    def __init__(self, hook: BoxHook):
        self.hook = hook

    def filter_files(self,
        files: XComArg,
        # A regex to match against any file name inside of a MatchedRecordFolder, resulting files will be linked.
        file_regex: re.Pattern = re.compile(r".*"),
        # This regex is matched first against the folder subpath within the matched record folder before linking all matching files (including none)
        folder_path_regex: re.Pattern = re.compile(r".*"),):
        """
        Creates a task that filters the Box files by regexp, and also applies a regexp against the
        folder path for each file.
        """

        @task
        def filter_files(files):
            result: List[BoxFile] = [f for f in files if folder_path_regex.match(f.folder_path) and file_regex.match(f.file_name)]

            return result

        return filter_files(files)

    def ingest_to_metis(self, project_name: str = "triage", bucket_name: str = "waiting_room") -> XComArg:
        pass

    def ingest_to_c4(
        self,
        # an XComArg that lists MatchedRecordFiles, usually the result of tail()
        listed_record_matches: XComArg,
        # By default, the ingestion process does a dry-run fetch for files.  set dry_run=False to actually download files to C4.
        dry_run=True,
    ) -> XComArg:
        """
        Creates a task that follows the listed set of record match folders and calls the given
        processor function to provide an UpdateRequest object as a result of matching files from
        the listed match record folder.

        The decorated function can receive arguments named metis, file, match, or template,
        each being either a MetisClient, File, MatchedRecordFolder, or Template respectively.

        eg:

        ```
        def my_processor(metis, template, match, file):
          with metis.open_file(file) as file_reader:
            csv_reader = csv.reader(file_reader)
            return match.as_update(csv_reader)

        matches = helpers.find_record_folders('rna_seq', rna_seq_folder_regex)
        listed_matches = helpers.list_match_folders(matches)
        helpers.process_and_link_matching_file(listed_matches, re.compile(r'.*gene_counts\.tsv$'), my_processor)
        ```
        """

        @task
        def do_link(listed_matches: List[Tuple[MatchedRecordFolder, List[File]]]):
            with self.hook.magma(
                project_name=get_project_name(), read_only=True
            ) as magma:
                models = magma.retrieve(get_project_name(), hide_templates=False).models
                with self.hook.metis(read_only=False) as metis:
                    for match, files in listed_matches:
                        for file in files:
                            if regex.match(file.file_name):
                                template: Template = None
                                if match.model_name in models:
                                    template = models[match.model_name].template

                                yield match, inject(
                                    processor,
                                    dict(
                                        metis=metis,
                                        file=file,
                                        match=match,
                                        template=template,
                                    ),
                                )

        return do_link(listed_record_matches)


def load_box_files_batch(
    box: Box, folder_name: str
) -> List[BoxFile]:
    return _load_box_files_batch(box, folder_name)


def _load_box_files_batch(
    box: Box,
    folder_name: str
) -> List[BoxFile]:
    context: Context = get_current_context()
    start, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for Box data from {start.isoformat(timespec='seconds')} to {end.isoformat(timespec='seconds')}"
    )

    files = box.tail(folder_name, start, end)

    return files