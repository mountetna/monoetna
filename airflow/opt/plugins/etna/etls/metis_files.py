import dataclasses
import re
from datetime import timezone
from datetime import timezone
from typing import Union, Literal, List, Optional, Tuple

import os.path
from airflow import AirflowException
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from serde import serialize, deserialize

from etna.etls.context import get_batch_range
from etna.hooks.etna import Metis, Folder, File

@serialize
@deserialize
@dataclasses.dataclass
class MatchedAtRoot:
    root_path: str
    match_file: Optional[File]
    match_folder: Optional[Folder]

    @property
    def root_base_name(self) -> str:
        return os.path.basename(self.root_path)

    @property
    def match_subpath(self) -> str:
        if self.match_file:
            return self.match_file.file_path[len(self.root_path) + 1:]
        if self.match_folder:
            return self.match_folder.folder_path[len(self.root_path) + 1:]
        return ''

    @property
    def match_full_path(self) -> str:
        return '/'.join(self.root_path, self.match_subpath)

def filter_by_root_directory(files_or_folders: Union[List[File], List[Folder]], directory_regex: re.Pattern) -> List[MatchedAtRoot]:
    result = []

    for file_or_folder in files_or_folders:
        path = ''
        file = None
        folder = None
        if isinstance(file_or_folder, File):
            path = file_or_folder.file_path
            file = file_or_folder
        if isinstance(file_or_folder, Folder):
            path = file_or_folder.folder_path
            folder = file_or_folder
        m = directory_regex.match(path)
        if not m:
            continue

        if m.pos != 0:
            continue

        # Don't catch matches that do not terminate on path segment
        if m.endpos < len(path):
            if path[m.endpos] != '/':
                continue

        result.append(
            MatchedAtRoot(
                root_path=path[:m.endpos],
                match_file=file,
                match_folder=folder,
            )
        )

    return result

# Increment this to force a reload of all metis cursor data.
metis_batch_loading_version = 1

def load_metis_files_and_folders_batch(
        metis: Metis,
        bucket_name: str,
        type: Union[Literal['file'], Literal['folder']],
        project_name: Optional[str]=None
    ) -> Union[List[File], List[Folder]]:
    context: Context = get_current_context()
    start, end = get_batch_range(context)

    project_name = project_name or metis.get_project_scope()
    if project_name is None:
        raise AirflowException("load_metis_files_and_folders_batch could not determine project_name from scope.")

    result = metis.find(
        project_name, bucket_name, [
            dict(
                type=type,
                attribute='updated_at',
                predicate='>=',
                value=start.replace(tzinfo=timezone.utc).isoformat(timespec='seconds')
            ),
            dict(
                type=type,
                attribute='updated_at',
                predicate='<=',
                value=end.replace(tzinfo=timezone.utc).isoformat(timespec='seconds')
            ),
        ]
    )

    if type == 'file':
        return result.files

    return result.folders

