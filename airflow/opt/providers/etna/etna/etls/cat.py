import logging
import os
import re
from typing import List, Optional

import paramiko

from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from airflow.models.xcom_arg import XComArg

from etna.hooks.etna import EtnaHook
from etna.etls.etl_task_batching import get_batch_range
from etna.hooks.cat import SftpEntry, Cat, CatHook
from etna.hooks.c4 import C4Hook
from etna.etls.remote_helpers_base import RemoteHelpersBase


class CatEtlHelpers(RemoteHelpersBase):
    hook: CatHook

    def __init__(self, hook: CatHook, magic_string: str):
        self.hook = hook
        self.log = logging.getLogger("airflow.task")
        self.magic_string = magic_string

        # paramiko.sftp_file.SFTPFile.MAX_REQUEST_SIZE = pow(2, 22) # 4MB per chunk

    def alert_slack(self,
        ingested_files: XComArg,
        ingested: bool,
        project_name: str,
        bucket_name: str,
        channel: str = "data-ingest-ping",
        target_system: str = "Metis",
        member_ids: Optional[List[str]] = None):
        """
        Sends a Slack message to the data-ingest-ping channel, notifying of
            the number of files uploaded.

        args:
            ingested_files: List of files
            ingested: bool, not really used, just helps control the flow of when messages are sent. Should be return value of helpers.
            project_name: str, project name for the message
            bucket_name: str, bucket name for the message
            channel: str, the Slack channel to post to, default data-ingest-ping
            target_system: str, where the files were ingested to, for the Slack message. Default of "Metis".
            member_ids: Optional[List[str]], list of Slack member ids to be notified of task completion.
        """

        @task
        def alert(ingested_files, ingested, project_name, bucket_name):
            self.alert(
                ingested_files,
                project_name,
                target_path=bucket_name,
                source_system="CAT",
                target_system=target_system,
                channel=channel,
                member_ids=member_ids)

        return alert(ingested_files, ingested, project_name, bucket_name)

    def ingest_to_c4(
        self,
        files: XComArg,
        folder_path: str = None,
        remove_magic_string: bool = True,
        batch_size: int = 5) -> XComArg:
        """
        Given a list of CAT files, will copy them to the given C4 path,
        mimicking the full directory structure from the CAT.

        args:
            files: List of files from tail_files or filter_files call
            folder_path: str, existing folder path to dump the files in. Default is C4 root_path configuration + folder structure on CAT.
            remove_magic_string: bool, remove the magic string from the file name. Default is True.
            batch_size: int, will save the cursor every X files that are ingested. Default is 5.
        """
        @task
        def ingest(files, folder_path):
            c4_hook = C4Hook.for_project()
            ingested_files = []
            with c4_hook.c4() as c4, self.hook.cat() as cat, c4.sftp() as c4_sftp:
                self.log.info(f"Attempting to upload {len(files)} files to C4")
                num_ingested = 0
                for file in files:
                    if cat.file_ingested_to_system("c4", file):
                        self.log.info(f"Skipping {file.name} because it has already been ingested.")
                        continue

                    final_file_name = file.name

                    if remove_magic_string:
                        final_file_name = file.name.replace(self.magic_string, "")

                    dest_path = os.path.join(folder_path or "", file.folder_path.replace(f"{cat._root_path()}/", ""))
                    self.log.info(f"Uploading {file.full_path} to {os.path.join(dest_path, final_file_name)}.")

                    c4.upload_file(
                        c4_sftp,
                        dest_path,
                        final_file_name,
                        file,
                        self.hook
                    )
                    cat.mark_file_as_ingested("c4", file)
                    self.log.info(f"Done ingesting {file.full_path}.")
                    num_ingested += 1
                    ingested_files.append(file)
                    if num_ingested % batch_size == 0:
                        cat.update_cursor("c4")
                cat.update_cursor("c4")
            return ingested_files

        return ingest(files, folder_path)

    def ingest_to_metis(
        self,
        files: XComArg,
        project_name: str = "triage",
        bucket_name: str = "waiting_room",
        folder_path: str = None,
        remove_magic_string: bool = True,
        batch_size: int = 5) -> XComArg:
        """
        Given a list of files, will copy them to the given Metis project_name and bucket_name,
        mimicking the full directory structure from the CAT.

        args:
            files: List of files from tail_files or filter_files call
            project_name: str, the target Metis project name. Default is `triage`
            bucket_name: str, the target Metis bucket name. Default is `waiting_room`
            folder_path: str, existing folder path to dump the files in. Default is Box hostname + folder structure in Box.
            remove_magic_string: bool, remove the magic string from the file name. Default is True.
            batch_size: int, will save the cursor every X files that are ingested. Default is 5.
        """
        @task
        def ingest(files, project_name, bucket_name, folder_path):
            etna_hook = EtnaHook.for_project(project_name)
            ingested_files = []
            with etna_hook.metis(project_name, read_only=False) as metis, self.hook.cat() as cat:
                self.log.info(f"Attempting to upload {len(files)} files to Metis")
                num_ingested = 0

                for file in files:
                    if cat.file_ingested_to_system("metis", file):
                        self.log.info(f"Skipping {file.name} because it has already been ingested.")
                        continue

                    with cat.retrieve_file(file) as file_handle:
                        dest_path = os.path.join(folder_path or "", file.folder_path.replace(f"{cat._root_path()}/", ""))

                        final_file_name = file.name

                        if remove_magic_string:
                            final_file_name = file.name.replace(self.magic_string, "")

                        self.log.info(f"Uploading {file.full_path} to {os.path.join(dest_path, final_file_name)}.")

                        self.handle_metis_ingest(
                            file_handle=file_handle,
                            folder_path=folder_path,
                            hostname=self.hook.connection.host,
                            flatten=False,
                            split_folder_name=None,
                            root_folder=cat._root_path(),
                            project_name=project_name,
                            bucket_name=bucket_name,
                            file=file,
                            metis=metis,
                            file_name_override=final_file_name
                        )
                    cat.mark_file_as_ingested("metis", file)
                    self.log.info(f"Done ingesting {file.full_path}.")
                    num_ingested += 1
                    ingested_files.append(file)
                    if num_ingested % batch_size == 0:
                        cat.update_cursor("metis")
                cat.update_cursor("metis")
            return ingested_files

        return ingest(files, project_name, bucket_name, folder_path)


def load_cat_files_batch(
    cat: Cat,
    magic_string: str,
    ignore_directories: List[str]
) -> List[SftpEntry]:
    return _load_cat_files_batch(
        cat,
        re.compile(f".*{magic_string}.*"),
        ignore_directories
    )


def _load_cat_files_batch(
    cat: Cat,
    magic_string: re.Pattern,
    ignore_directories: List[str]
) -> List[SftpEntry]:
    context: Context = get_current_context()
    _, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for CAT data that has not been ingested as of {end.isoformat(timespec='seconds')}"
    )

    files = cat.tail(
        magic_string=magic_string,
        ignore_directories=ignore_directories
    )

    return files
