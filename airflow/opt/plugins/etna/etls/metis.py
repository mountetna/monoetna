import dataclasses
import re
from datetime import timezone
from datetime import timezone
from typing import Union, Literal, List, Optional, Tuple, Callable

import os.path
from airflow.exceptions import AirflowException
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context, task
from serde import serialize, deserialize

from etna.etls.context import get_batch_range
from etna.hooks.etna import Metis, Folder, File, Magma


@serialize
@deserialize
@dataclasses.dataclass
class MatchedAtRoot:
    root_path: str
    match_file: Optional[File]
    match_folder: Optional[Folder]

    @property
    def updated_at(self):
        if self.match_file:
            return self.match_file.updated_at
        return self.match_folder.updated_at

    @property
    def project_name(self) -> str:
        if self.match_file:
            return self.match_file.project_name
        return self.match_folder.project_name

    @property
    def bucket_name(self) -> str:
        if self.match_file:
            return self.match_file.bucket_name
        return self.match_folder.bucket_name

    @property
    def folder_path(self) -> str:
        if self.match_file:
            return os.path.dirname(self.match_file.file_path)
        return self.match_folder.folder_path

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
        return '/'.join([self.root_path, self.match_subpath])

def link(model_name, attribute_name, dry_run=True):
    pass

def list_contents_of_matches(metis: Metis, matching: List[MatchedAtRoot]) -> List[Tuple[MatchedAtRoot, File]]:
    return [
        (m, file)
        for file in metis.list_folder(m.project_name, m.bucket_name, m.folder_path).files
        for m in matching
    ]


def filter_by_exists_in_timur(
        magma: Magma,
        matched: List[MatchedAtRoot],
        model_name: str,
        corrected_record_name: Callable[[str], str] = lambda x: x
):
    if not matched:
        return []

    project_name = matched[0].project_name
    record_names = [corrected_record_name(m.root_base_name) for m in matched]
    response = magma.retrieve(project_name,
                              model_name=model_name,
                              attribute_names='identifier',
                              record_names=record_names)
    return [m for m in matched if corrected_record_name(m.root_base_name) in response.models[model_name].documents]


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

        start, end = m.span()

        if start != 0:
            continue

        print(m, start, end)

        # Don't catch matches that do not terminate on path segment
        if end < len(path):
            if path[end] != '/':
                continue

        result.append(
            MatchedAtRoot(
                root_path=path[:end],
                match_file=file,
                match_folder=folder,
            )
        )

    return result

# Increment this to force a reload of all metis cursor data.
metis_batch_loading_version = 1

def load_metis_files_batch(metis: Metis, bucket_name: str, project_name: Optional[str]=None) -> List[File]:
    return _load_metis_files_and_folders_batch(metis, bucket_name, 'file', project_name)

def load_metis_folders_batch(metis: Metis, bucket_name: str, project_name: Optional[str]=None) -> List[Folder]:
    return _load_metis_files_and_folders_batch(metis, bucket_name, 'folder', project_name)

def _load_metis_files_and_folders_batch(
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

