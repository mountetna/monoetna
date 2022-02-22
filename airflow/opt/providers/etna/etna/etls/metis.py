import dataclasses
import functools
import io
import logging
import os.path
import re
from logging import Logger
from typing import Union, Literal, List, Optional, Tuple, Callable, Dict, Any, Iterable

from airflow.decorators import task
from airflow.exceptions import AirflowException
from airflow.models.taskinstance import Context
from airflow.models.xcom_arg import XComArg
from airflow.operators.python import get_current_context
from serde import serialize, deserialize
from serde.json import from_json

from etna.dags.project_name import get_project_name
from etna.etls.etl_task_batching import get_batch_range
from etna.hooks.etna import File, Folder, UpdateRequest, EtnaHook, Magma, Metis, Model
from etna.operators import DockerOperatorBase
from etna.utils.iterables import batch_iterable
from etna.xcom.etna_xcom import pickled


@serialize
@deserialize
@dataclasses.dataclass
class MatchedRecordFolder:
    root_path: str
    record_name: str
    model_name: str
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

    def as_update(self, attrs: Dict[str, Any]):
        req = UpdateRequest(project_name=self.project_name)
        req.update_record(self.model_name, self.record_name, attrs)
        return req

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
    def match_folder_subpath(self) -> str:
        return self.folder_path[len(self.root_path) + 1 :]

    @property
    def match_subpath(self) -> str:
        if self.match_file:
            return self.match_file.file_path[len(self.root_path) + 1 :]
        if self.match_folder:
            return self.match_folder.folder_path[len(self.root_path) + 1 :]
        return ""

    @property
    def match_full_path(self) -> str:
        if not self.match_subpath:
            return self.root_path
        return "/".join([self.root_path, self.match_subpath])


class link:
    task_id: Optional[str]
    log: Logger
    version: int
    dry_run: bool
    _hook: Optional[EtnaHook]
    attribute_name: Optional[str]

    def __init__(
        self,
        attribute_name: Optional[str] = None,
        dry_run=True,
        hook: Optional[EtnaHook] = None,
        task_id: Optional[str] = None,
    ):
        self.task_id = task_id
        self._hook = hook
        self.dry_run = dry_run
        self.attribute_name = attribute_name
        self.log = logging.getLogger("airflow.task")

    @property
    def project_name(self) -> str:
        return get_project_name()

    @property
    def hook(self):
        return self._hook or EtnaHook.for_project(self.project_name)

    def __call__(self, fn: Callable):
        @functools.wraps(fn)
        def new_task(*args, **kwds):
            result: List[MatchedRecordFolder] = []
            with self.hook.magma(read_only=False) as magma:
                self.log.info("retrieving model...")
                response = magma.retrieve(
                    project_name=self.project_name,
                    model_name="all",
                    attribute_names=[self.attribute_name]
                    if self.attribute_name
                    else "all",
                    record_names=[],
                    hide_templates=False,
                )

                self.log.info("running linking function...")
                for cur_batch in batch_iterable(fn(*args, **kwds), 50):
                    result.extend(
                        self._process_link_batch(magma, response.models, cur_batch)
                    )

            return pickled(result)

        if self.task_id:
            fn.__name__ = self.task_id

        new_task.__name__ = f"{fn.__name__}_{'_dry_run' if self.dry_run else ''}"

        return task()(new_task)

    def _process_link_batch(
        self,
        magma: Magma,
        models: Dict[str, Model],
        batch: Iterable[Tuple[MatchedRecordFolder, Any]],
    ) -> List[MatchedRecordFolder]:
        result: List[MatchedRecordFolder] = []
        update: UpdateRequest = UpdateRequest(
            revisions={}, project_name=self.project_name, dry_run=self.dry_run
        )
        for match, value in batch:
            if self.dry_run:
                pass

            if isinstance(value, UpdateRequest):
                pass
            elif self.attribute_name is not None:
                value = match.as_update({self.attribute_name: value})
            else:
                raise AirflowException(
                    f"Cannot link {value}, attribute_name must be set or UpdateRequest returned from linking function."
                )

            for error_message in value.validate(models):
                raise AirflowException(error_message)

            update.extend(value)
            result.append(match)

        self.log.info("Executing batched magma update")
        magma.update(update)
        return result


def list_contents_of_matches(
    metis: Metis, matching: List[MatchedRecordFolder]
) -> List[Tuple[MatchedRecordFolder, List[File]]]:
    return [
        (m, metis.list_folder(m.project_name, m.bucket_name, m.folder_path).files)
        for m in matching
    ]


class MetisEtlHelpers:
    hook: EtnaHook

    def __init__(self, tail_folders: XComArg, tail_files: XComArg, hook: EtnaHook):
        self.hook = hook
        self.tail_files = tail_files
        self.tail_folders = tail_folders

    def prepare_task_token(self) -> Callable[[bool], XComArg]:
        @task
        def prepare_task_token(read_only):
            return self.hook.get_task_auth(read_only=read_only).token

        return prepare_task_token

    def process_and_link_matching_file(
        self,
        listed_record_matches: XComArg,
        regex: re.Pattern,
        processor: Callable[[Metis, MatchedRecordFolder, File], UpdateRequest],
        dry_run=True,
    ) -> XComArg:
        @link(dry_run=dry_run, task_id=processor.__name__)
        def do_link(listed_matches: List[Tuple[MatchedRecordFolder, List[File]]]):
            with self.hook.metis(read_only=False) as metis:
                for match, files in listed_matches:
                    for file in files:
                        if regex.match(file.file_name):
                            yield match, processor(metis, match, file)

        return do_link(listed_record_matches)

    def link_matching_file(
        self,
        listed_record_matches: XComArg,
        attr_name: str,
        regex: re.Pattern,
        dry_run=True,
    ) -> XComArg:
        @link(attr_name, dry_run=dry_run, task_id=f"link_{attr_name}")
        def do_link(listed_matches: List[Tuple[MatchedRecordFolder, List[File]]]):
            for match, files in listed_matches:
                for file in files:
                    if regex.match(file.file_name):
                        yield match, file

        return do_link(listed_record_matches)

    def link_matching_files(
        self,
        listed_record_matches: XComArg,
        attr_name: str,
        file_regex: re.Pattern,
        folder_path_regex: re.Pattern = re.compile(r"^$"),
        dry_run=True,
    ) -> XComArg:
        @link(attr_name, dry_run=dry_run, task_id=f"link_{attr_name}")
        def do_link(listed_matches: List[Tuple[MatchedRecordFolder, List[File]]]):
            for match, files in listed_matches:
                if folder_path_regex.match(match.match_folder_subpath):
                    matched_files = [f for f in files if file_regex.match(f.file_name)]
                    yield match, matched_files

        return do_link(listed_record_matches)

    def find_record_folders(
        self,
        model_name: str,
        regex: re.Pattern,
        corrected_record_name: Callable[[str], str] = lambda x: x,
    ) -> XComArg:
        @task
        def find_record_folders(folders, files):
            return pickled(
                filter_by_record_directory(
                    folders + files,
                    regex,
                    model_name=model_name,
                    corrected_record_name=corrected_record_name,
                )
            )

        return find_record_folders(self.tail_folders, self.tail_files)

    def filter_by_timur(self, matches: XComArg) -> XComArg:
        @task
        def filter_by_timur(matches):
            with self.hook.magma() as magma:
                return pickled(filter_by_exists_in_timur(magma, matches))

        return filter_by_timur(matches)

    def list_match_folders(self, matches: XComArg):
        @task
        def list_match_folders(matches):
            with self.hook.metis() as metis:
                return pickled(list_contents_of_matches(metis, matches))

        return list_match_folders(matches)


def filter_by_exists_in_timur(
    magma: Magma,
    matched: List[MatchedRecordFolder],
):
    if not matched:
        return []

    project_name = matched[0].project_name
    matched_by_model_name = {}
    for m in matched:
        matched_by_model_name.setdefault(m.model_name, []).append(m)

    return [
        m
        for model_name, matches in matched_by_model_name
        for response in [
            magma.retrieve(
                project_name,
                model_name=model_name,
                attribute_names="identifier",
                record_names=[m.record_name for m in matches],
            )
        ]
        for m in matches
        if m.model_name in response.models[model_name].documents
    ]


def filter_by_record_directory(
    files_or_folders: List[Union[File, Folder]],
    directory_regex: re.Pattern,
    model_name: str,
    corrected_record_name: Callable[[str], str] = lambda x: x,
) -> List[MatchedRecordFolder]:
    result: Dict[str, MatchedRecordFolder] = {}

    for file_or_folder in files_or_folders:
        path = ""
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
            if path[end] != "/":
                continue

        root_path = path[:end]
        match = MatchedRecordFolder(
            root_path=root_path,
            record_name=corrected_record_name(os.path.basename(root_path)),
            match_file=file,
            match_folder=folder,
            model_name=model_name,
        )

        if match.folder_path in result:
            continue

        result[match.folder_path] = match

    return list(result.values())


# Increment this to force a reload of all metis cursor data.
metis_batch_loading_version = 1


def load_metis_files_batch(
    metis: Metis, bucket_name: str, project_name: Optional[str] = None
) -> List[File]:
    return _load_metis_files_and_folders_batch(metis, bucket_name, "file", project_name)


def load_metis_folders_batch(
    metis: Metis, bucket_name: str, project_name: Optional[str] = None
) -> List[Folder]:
    return _load_metis_files_and_folders_batch(
        metis, bucket_name, "folder", project_name
    )


def _load_metis_files_and_folders_batch(
    metis: Metis,
    bucket_name: str,
    type: Union[Literal["file"], Literal["folder"]],
    project_name: Optional[str] = None,
) -> Union[List[File], List[Folder]]:
    context: Context = get_current_context()
    start, end = get_batch_range(context)

    project_name = project_name or metis.get_project_scope()
    if project_name is None:
        raise AirflowException(
            "load_metis_files_and_folders_batch could not determine project_name from scope."
        )

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for metis data from {start.isoformat(timespec='seconds')} to {end.isoformat(timespec='seconds')}"
    )

    response = metis.find(
        project_name,
        bucket_name,
        [
            dict(
                type=type,
                attribute="updated_at",
                predicate=">=",
                value=start.isoformat(timespec="seconds"),
            ),
            dict(
                type=type,
                attribute="updated_at",
                predicate="<=",
                value=end.isoformat(timespec="seconds"),
            ),
        ],
    )

    if type == "file":
        return response.files
    else:
        return response.folders
