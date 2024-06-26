from typing import List, Optional
import os
import logging
import socket
import time

from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from airflow.models.xcom_arg import XComArg

from etna.etls.etl_task_batching import get_batch_range
from etna.etls.remote_helpers_base import RemoteHelpersBase
from etna.hooks.box import BoxHook, FtpEntry, Box
from etna.hooks.etna import EtnaHook


class BoxEtlHelpers(RemoteHelpersBase):
    hook: BoxHook

    def __init__(self, hook: BoxHook, box_folder: str):
        self.hook = hook
        self.log = logging.getLogger("airflow.task")
        self.box_folder = box_folder

    def alert_slack(self,
        files: XComArg,
        ingested: bool,
        project_name: str,
        bucket_name: str,
        channel: str = "data-ingest-ping",
        member_ids: Optional[List[str]] = None,
        retries: Optional[int] = 10,
        depends_on_past: Optional[bool] = False):
        """
        Sends a Slack message to the data-ingest-ping channel, notifying of
            the number of files uploaded.

        args:
            files: List of files
            ingested: bool, not really used, just helps control the flow of when messages are sent. Should be return value of helpers.ingest_to_metis.
            project_name: str, project name for the message
            bucket_name: str, bucket name for the message
            channel: str, the Slack channel to post to, default data-ingest-ping
            member_ids: Optional[List[str]], list of Slack member ids to be notified of task completion.
            retries: Optional[int], number of retries for this task. Default of 10.
            depends_on_past: Apache Airflow depends_on_past flag -- can only run if the previous run succeeded. Default of False.
        """

        @task(retries=retries, depends_on_past=depends_on_past)
        def alert(files, ingested, project_name, bucket_name):
            self.alert(
                files,
                project_name,
                target_path=bucket_name,
                source_system="Box",
                target_system="Metis",
                channel=channel,
                member_ids=member_ids)

        return alert(files, ingested, project_name, bucket_name)

    def ingest_to_metis(
        self,
        files: XComArg,
        project_name: str = "triage",
        bucket_name: str = "waiting_room",
        folder_path: str = None,
        flatten: bool = True,
        clean_up: bool = False,
        split_folder_name: str = None,
        batch_size: int = 5,
        retries: Optional[int] = 10,
        depends_on_past: Optional[bool] = False) -> XComArg:
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
            split_folder_name: str, if flatten=False, the folder name after which to copy the structure from Box. Generally would match what you use in filter_files() for folder_path_regex.
            batch_size: int, number of files to save in the cursor, at a time. Default is 5.
            retries: Optional[int], number of retries for this task. Default of 10.
            depends_on_past: Apache Airflow depends_on_past flag -- can only run if the previous run succeeded. Default of False.
        """
        @task(retries=retries, depends_on_past=depends_on_past)
        def ingest(files, project_name, bucket_name, folder_path, batch_size):
            etna_hook = EtnaHook.for_project(project_name)
            with etna_hook.metis(project_name, read_only=False) as metis, self.hook.box() as box:
                self.log.info(f"Attempting to upload {len(files)} files to Metis")
                num_ingested = 0
                for file in files:
                    if box.file_ingested_to_system(file):
                        self.log.info(f"Skipping {file.name} because it has already been ingested.")
                        continue

                    with box.ftps() as ftps:
                        sock = box.retrieve_file(ftps, file)

                        with sock.makefile(mode='rb') as connection:
                            self.handle_metis_ingest(
                                file_handle=connection,
                                folder_path=folder_path,
                                hostname=self.hook.connection.host,
                                flatten=flatten,
                                split_folder_name=split_folder_name,
                                root_folder=self.box_folder,
                                project_name=project_name,
                                bucket_name=bucket_name,
                                file=file,
                                metis=metis
                            )

                            if clean_up:
                                box.remove_file(ftps, file)
                                self.log.info("Removed the original file from Box.")
                            box.mark_file_as_ingested(file)
                            self.log.info(f"Done ingesting {file.full_path}.")
                            num_ingested += 1
                            if num_ingested % batch_size == 0:
                                box.update_cursor()
                        sock.shutdown(socket.SHUT_RDWR)
                        sock.close()
                box.update_cursor()

        return ingest(files, project_name, bucket_name, folder_path, batch_size)


def load_box_files_batch(
    box: Box, folder_name: str
) -> List[FtpEntry]:
    return _load_box_files_batch(box, folder_name)


def _load_box_files_batch(
    box: Box,
    folder_name: str
) -> List[FtpEntry]:
    context: Context = get_current_context()
    _, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for Box data that has not been ingested as of {end.isoformat(timespec='seconds')}"
    )

    files = box.tail(folder_name)

    return files
