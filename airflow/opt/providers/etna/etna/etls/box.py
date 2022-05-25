from typing import List, Optional
import os
import re
import logging
import socket
import time

from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from airflow.models.xcom_arg import XComArg
from airflow.providers.slack.operators.slack_webhook import SlackWebhookOperator

from etna.etls.etl_task_batching import get_batch_range

from etna.hooks.box import BoxHook, FtpEntry, Box
from etna.hooks.etna import EtnaHook


class BoxEtlHelpers:
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
        member_ids: Optional[List[str]] = None):
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
        """

        @task
        def alert(files, ingested, project_name, bucket_name):
            if len(files) > 0:
                msg = "\n".join([f"Finished uploading {len(files)} files from Box to Metis for {project_name}. Please check the {bucket_name} bucket."] + [f.full_path for f in files])

                SlackWebhookOperator(
                    task_id=f"notify_slack_{project_name}_{bucket_name}_box_ingest",
                    username="Airflow",
                    channel=channel,
                    http_conn_id='slack-api',
                    message=msg
                ).execute(context=None)

                if member_ids is not None and len(member_ids) > 0:
                    user_mentions = [f"<@{m_id}>" for m_id in member_ids]
                    notify_msg = f"{' '.join(user_mentions)} :point_up_2:"
                    SlackWebhookOperator(
                        task_id=f"notify_slack_users_{project_name}_{bucket_name}_box_ingest",
                        username="Airflow",
                        channel=channel,
                        http_conn_id='slack-api',
                        message=notify_msg
                    ).execute(context=None)

        return alert(files, ingested, project_name, bucket_name)

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
        clean_up: bool = False,
        split_folder_name: str = None) -> XComArg:
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
        """
        @task
        def ingest(files, project_name, bucket_name, folder_path):
            etna_hook = EtnaHook.for_project(project_name)
            with etna_hook.metis(project_name, read_only=False) as metis, self.hook.box() as box:
                self.log.info(f"Attempting to upload {len(files)} files to Metis")
                for file in files:
                    with box.ftps() as ftps:
                        sock = box.retrieve_file(ftps, file)

                        with sock.makefile(mode='rb') as connection:
                            parts = []
                            if folder_path is None:
                                parts.append(self.hook.connection.host)
                            else:
                                parts.append(folder_path)

                            if flatten:
                                parts.append(file.name)
                            elif split_folder_name is not None:
                                rel_path = file.rel_path
                                parts.append(rel_path.split(f"/{split_folder_name}/")[-1])
                            else:
                                rel_path = file.rel_path
                                parts.append(rel_path.split(f"/{self.box_folder}/")[-1])

                            dest_path = os.path.join(*parts)

                            self.log.info(f"Uploading {file.full_path} to {dest_path}.")

                            should_log = True
                            for blob in metis.upload_file(
                                project_name,
                                bucket_name,
                                dest_path,
                                connection,
                                file.size
                            ):
                                # Only log every 5 seconds, to save log space...
                                time_check = int(time.time())
                                if time_check % 5 == 0 and should_log:
                                    self.log.info("Uploading blob...")
                                    should_log = False
                                elif time_check % 5 != 0 and not should_log:
                                    should_log = True

                            if clean_up:
                                box.remove_file(ftps, file)
                                self.log.info("Removed the original file from Box.")
                            box.mark_file_as_ingested(file)
                            self.log.info(f"Done ingesting {file.full_path}.")
                        sock.shutdown(socket.SHUT_RDWR)
                        sock.close()
                box.update_cursor()

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
    _, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for Box data that has not been ingested as of {end.isoformat(timespec='seconds')}"
    )

    files = box.tail(folder_name)

    return files
