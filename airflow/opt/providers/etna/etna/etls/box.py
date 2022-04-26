from typing import List
import os
import re
import logging
import time

from airflow.operators.python import get_current_context

from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from airflow.models.xcom_arg import XComArg

from etna.etls.etl_task_batching import get_batch_range

from etna.hooks.box import BoxHook, FtpEntry, Box
from etna.hooks.etna import EtnaHook


class BoxEtlHelpers:
    hook: BoxHook

    def __init__(self, hook: BoxHook):
        self.hook = hook
        self.log = logging.getLogger("airflow.task")

    def filter_files(self,
        files: XComArg,
        # A regex to match against any file name, resulting files will be linked.
        file_name_regex: re.Pattern = re.compile(r".*"),
        # This regex is matched against the folder subpath of the file
        folder_path_regex: re.Pattern = re.compile(r".*"),):
        """
        Creates a task that filters the Box file names by regexp, and can also apply a regexp against the
        folder path for each file for further filtering.

        args:
            files: List of files from tail_files or previous filter_files call
            file_name_regex: re.Pattern, i.e. re.compile(r".*\.fcs")
            folder_path_regex: re.Pattern, i.e. re.compile(r".*ipi.*")
        """

        @task
        def filter_files(files):
            result: List[FtpEntry] = [f for f in files if folder_path_regex.match(f.folder_path) and file_name_regex.match(f.name)]

            return result

        return filter_files(files)

    def ingest_to_metis(
        self,
        files: XComArg,
        project_name: str = "triage",
        bucket_name: str = "waiting_room",
        folder_path: str = None,
        flatten: bool = True,
        clean_up: bool = False) -> XComArg:
        """
        Given a list of Box files, will copy them to the given Metis project_name and bucket_name,
        mimicking the full directory structure from Box.

        args:
            files: List of files from tail_files or filter_files call
            project_name: str, the target Metis project name. Default is `triage`
            bucket_name: str, the target Metis bucket name. Default is `waiting_room`
            folder_path: str, existing folder path to dump the files in. Default is Box hostname + folder structure in Box.
            flatten: bool, to flatten the Box folder structure or maintain it. Default is True.
            clean_up: bool, to remove the file from Box after ingest. Default is False.
        """
        @task
        def ingest(files, project_name, bucket_name, folder_path):
            etna_hook = EtnaHook.for_project(project_name)
            with etna_hook.metis(project_name, read_only=False) as metis, self.hook.box() as box, box.ftps() as ftps:
                self.log.info(f"Attempting to upload {len(files)} files to Metis")
                for file in files:
                    sock, size = box.retrieve_file(ftps, file)

                    with sock as connection:
                        parts = []
                        if folder_path is None:
                            parts.append(self.hook.connection.host)
                        else:
                            parts.append(folder_path)

                        if flatten:
                            parts.append(file.name)
                        else:
                            parts.append(file.rel_path)

                        dest_path = os.path.join(*parts)

                        self.log.info(f"Uploading {file.full_path} to {dest_path}.")

                        for blob in metis.upload_file(
                            project_name,
                            bucket_name,
                            dest_path,
                            connection,
                            size or file.size
                        ):
                            # Only log every 5 seconds, to save log space...
                            if int(time.time()) % 5 == 0:
                                self.log.info("Uploading blob...")
                        if clean_up:
                            box.remove_file(ftps, file)
                        self.log.info(f"Done ingesting {file.full_path}.")

        return ingest(files, project_name, bucket_name, folder_path)


def load_box_files_batch(
    box: Box, folder_name: str
) -> List[FtpEntry]:
    return _load_box_files_batch(box, folder_name)


def _load_box_files_batch(
    box: Box,
    folder_name: str
) -> List[FtpEntry]:
    context: Context = get_current_context()
    start, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for Box data from {start.isoformat(timespec='seconds')} to {end.isoformat(timespec='seconds')}"
    )

    files = box.tail(folder_name, start, end)

    return files
