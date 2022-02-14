import dataclasses
import functools
import itertools
import logging
import re
from datetime import timezone
from typing import Union, Literal, List, Optional, Tuple, Callable, Dict

import os.path

from airflow.decorators import task
from airflow.exceptions import AirflowException
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context
from serde import serialize, deserialize

from etna.dags.project_name import get_project_name
from etna.etls.etl_task_batching import get_batch_range
from etna.hooks.etna import File, Folder, UpdateRequest, EtnaHook, Magma, Metis
from etna.utils.iterables import batch_iterable


@serialize
@deserialize
@dataclasses.dataclass
class MatchedAtRoot:
    root_path: str
    record_name: str
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
    def match_subpath(self) -> str:
        if self.match_file:
            return self.match_file.file_path[len(self.root_path) + 1:]
        if self.match_folder:
            return self.match_folder.folder_path[len(self.root_path) + 1:]
        return ''

    @property
    def match_full_path(self) -> str:
        if not self.match_subpath:
            return self.root_path
        return '/'.join([self.root_path, self.match_subpath])

def link(model_name, attribute_name, dry_run=True, hook: Optional[EtnaHook] = None):
    def wrapper(fn: Callable):
        @functools.wraps(fn)
        def new_task(*args, **kwds):
            nonlocal hook
            project_name = get_project_name()

            log = logging.getLogger('airflow.task')

            if hook is None:
                hook = EtnaHook.for_project(project_name)

            with hook.magma(read_only=False) as magma:
                log.info('retrieving model...')
                response = magma.retrieve(
                    project_name=project_name, model_name=model_name,
                    attribute_names=[attribute_name], record_names=[],
                    hide_templates=False
                )

                template = response.models[model_name].template

                if attribute_name not in template.attributes:
                    raise AirflowException(f"Attribute '{attribute_name}' could not be linked in '{model_name}', as it does not exist.")
                attribute = template.attributes[attribute_name]

                log.info('running linking function...')
                for cur_batch in batch_iterable(fn(*args, **kwds), 50):
                    update: UpdateRequest = UpdateRequest(revisions={}, project_name=project_name, dry_run=dry_run)
                    for match, files in cur_batch:
                        if match.record_name in update.revisions:
                            log.warning(f"Found multiple matches for record '{match.record_name}': {' '.join(f.file_path for f in files)}")

                        record = update.update_record(model_name, match.record_name)

                        if attribute.attribute_type == "file_collection":
                            record.setdefault(attribute_name, [])

                        for file in files:
                            log.info(f"Found match for record #{match.record_name} on {file.file_path}")
                            if attribute.attribute_type == "file":
                                if attribute_name in record:
                                    log.warning("Multiple files for single file collection!")
                                record[attribute_name] = file.as_magma_file_attribute
                            elif attribute.attribute_type == "file_collection":
                                record[attribute_name].append(file.as_magma_file_attribute)
                        update.update_record(model_name, match.record_name)

                    log.info("Executing magma update")
                    magma.update(update)

        new_task.__name__ = f"{fn.__name__}{'_dry_run' if dry_run else ''}"
        return task(do_xcom_push=False)(new_task)

    return wrapper

def list_contents_of_matches(metis: Metis, matching: List[MatchedAtRoot]) -> List[Tuple[MatchedAtRoot, List[File]]]:
    return [
        (m, metis.list_folder(m.project_name, m.bucket_name, m.folder_path).files)
        for m in matching
    ]


def filter_by_exists_in_timur(
        magma: Magma,
        matched: List[MatchedAtRoot],
        model_name: str,
):
    if not matched:
        return []

    project_name = matched[0].project_name
    record_names = [m.record_name for m in matched]

    response = magma.retrieve(project_name,
                              model_name=model_name,
                              attribute_names='identifier',
                              record_names=record_names)
    return [m for m in matched if m.record_name in response.models[model_name].documents]


def filter_by_record_directory(
        files_or_folders: List[Union[File, Folder]],
        directory_regex: re.Pattern,
        corrected_record_name: Callable[[str], str] = lambda x: x,
) -> List[MatchedAtRoot]:
    result: Dict[str, MatchedAtRoot] = {}

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

        # Don't catch matches that do not terminate on path segment
        if end < len(path):
            if path[end] != '/':
                continue

        root_path = path[:end]
        if root_path in result:
            continue

        result[root_path] = MatchedAtRoot(
            root_path=root_path,
            record_name=corrected_record_name(os.path.basename(root_path)),
            match_file=file,
            match_folder=folder,
        )

    return list(result.values())

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

    response = metis.find(
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
        return response.files
    else:
        return response.folders
