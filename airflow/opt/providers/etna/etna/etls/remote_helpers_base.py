# Helper class to send Slack messages

import io
import os
import re
import time
from typing import List, Optional

from airflow.decorators import task
from airflow.models.xcom_arg import XComArg
from airflow.providers.slack.operators.slack_webhook import SlackWebhookOperator

from etna.hooks.hook_helpers import RemoteFileBase
from etna.hooks.ssh_base import SSHBase

from mountetna import Metis
class RemoteHelpersBase:
    def filter_files(self,
        files: XComArg,
        # A regex to match against any file name, resulting files will be linked.
        file_name_regex: re.Pattern = re.compile(r".*"),
        # This regex is matched against the folder subpath of the file
        folder_path_regex: re.Pattern = re.compile(r".*"),):
        """
        Creates a task that filters the given file names by regexp, and can also apply a regexp against the
        folder path for each file for further filtering.

        args:
            files: List of files from tail_files or previous filter_files call
            file_name_regex: re.Pattern, i.e. re.compile(r".*\.fcs")
            folder_path_regex: re.Pattern, i.e. re.compile(r".*ipi.*")
        """

        @task
        def filter_files(files):
            result: List[RemoteFileBase] = [f for f in files if folder_path_regex.match(f.folder_path) and file_name_regex.match(f.name)]

            return result

        return filter_files(files)

    def alert(self,
        files: XComArg,
        project_name: str,
        target_path: str,
        source_system: str,
        target_system: str,
        channel: str = "data-ingest-ping",
        member_ids: Optional[List[str]] = None):
        """
        Sends a Slack message to the data-ingest-ping channel, notifying of
            the number of files uploaded.

        args:
            files: List of files
            project_name: str, project name for the message
            target_path: str, target_path, i.e. the bucket or directory path where files were copied to
            source_system: str, where the files are coming from
            target_system: str, where the files were copied to
            channel: str, the Slack channel to post to, default data-ingest-ping
            member_ids: Optional[List[str]], list of Slack member ids to be notified of task completion.
        """
        if len(files) > 0:
            msg = "\n".join([f"Finished uploading {len(files)} files from {source_system} to {target_system} for {project_name}. Please check {target_path} bucket / path."] + [f.full_path for f in files])

            SlackWebhookOperator(
                task_id=f"notify_slack_{project_name}_{target_path}_{source_system}_ingest",
                username="Airflow",
                channel=channel,
                http_conn_id='slack-api',
                message=msg
            ).execute(context=None)

            if member_ids is not None and len(member_ids) > 0:
                user_mentions = [f"<@{m_id}>" for m_id in member_ids]
                notify_msg = f"{' '.join(user_mentions)} :point_up_2:"
                SlackWebhookOperator(
                    task_id=f"notify_slack_users_{project_name}_{target_path}_{source_system}_ingest",
                    username="Airflow",
                    channel=channel,
                    http_conn_id='slack-api',
                    message=notify_msg
                ).execute(context=None)

    def handle_metis_ingest(
        self,
        file_handle: io.BytesIO,
        folder_path: str,
        hostname: str,
        flatten: bool,
        split_folder_name: str,
        root_folder: str,
        project_name: str,
        bucket_name: str,
        file: RemoteFileBase,
        metis: Metis):
        parts = []
        if folder_path is None:
            parts.append(hostname)
        else:
            parts.append(folder_path)

        if flatten:
            parts.append(file.name)
        elif split_folder_name is not None:
            rel_path = file.rel_path
            parts.append(rel_path.split(f"/{split_folder_name}/")[-1])
        else:
            rel_path = file.rel_path
            parts.append(rel_path.split(f"{root_folder}/")[-1])

        dest_path = os.path.join(*parts)

        self.log.info(f"Uploading {file.full_path} to {dest_path}.")

        should_log = True
        for blob in metis.upload_file(
            project_name,
            bucket_name,
            dest_path,
            file_handle,
            file.size
        ):
            # Only log every 5 seconds, to save log space...
            time_check = int(time.time())
            if time_check % 5 == 0 and should_log:
                self.log.info("Uploading blob...")
                should_log = False
            elif time_check % 5 != 0 and not should_log:
                should_log = True

