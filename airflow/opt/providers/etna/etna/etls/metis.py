import dataclasses
import difflib
import functools
import io
import itertools
import json
import logging
import os.path
import re
from logging import Logger
from typing import Union, Literal, List, Optional, Tuple, Callable, Dict, Any, Iterable

from airflow import DAG
from airflow.decorators import task
from airflow.exceptions import AirflowException
from airflow.models.taskinstance import Context, TaskInstance
from airflow.models.xcom_arg import XComArg
from airflow.operators.python import get_current_context
from serde import serialize, deserialize
from serde.json import from_json

from etna.dags.project_name import get_project_name
from etna.etls.etl_task_batching import get_batch_range
from etna.hooks.etna import (
    File,
    Folder,
    UpdateRequest,
    EtnaHook,
    Magma,
    Metis,
    Model,
    Template, RetrievalResponse,
)
from etna.operators import DockerOperatorBase
from etna.utils.inject import inject
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
    def folder_id(self) -> int:
        """
        The folder id that was matched inside of a record folder.
        Note, this is NOT necessarily the id of the record folder itself, but either
        1.  the id of the folder found inside the record folder match or
        2.  the id of the parent folder of the file found inside the record folder (which, could be the record folder)
        """
        if self.match_file:
            return self.match_file.folder_id
        return self.match_folder.id

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
    debug_bucket: Optional[str]

    def __init__(
        self,
        attribute_name: Optional[str] = None,
        dry_run=True,
        hook: Optional[EtnaHook] = None,
        task_id: Optional[str] = None,
        debug_bucket: Optional[str] = None,
    ):
        self.task_id = task_id
        self._hook = hook
        self.dry_run = dry_run
        self.attribute_name = attribute_name
        self.log = logging.getLogger("airflow.task")
        self.debug_bucket = debug_bucket

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

        self._upload_debug(update, models)
        self.log.info("Executing batched magma update")
        magma.update(update)
        return result

    def _upload_debug(self, update: UpdateRequest, models: Dict[str, Model]):
        context = get_current_context()
        dag: DAG = context['dag']
        ti: TaskInstance = context['ti']

        if not self.debug_bucket:
            return

        self.log.info('Starting debug upload')
        records: RetrievalResponse = RetrievalResponse()
        with self.hook.magma() as magma:
            with self.hook.metis(read_only=False) as metis:
                for model_name, revisions in update.revisions.items():
                    records.extend(magma.retrieve(get_project_name(), model_name=model_name, record_names=list(revisions.keys())))
                self.log.info('prepared all revisions by model_name')

                for identifier, revision in revisions.items():
                    existing = records.models[model_name].documents.get(identifier, {})
                    # Compare the subset of attributes that actually exist in this revision
                    # against what might already exist for that record.
                    # Order the attributes based on the template to ensure cleaner diff as well.
                    exisiting_diffable = dict()
                    revision_diffable = dict()

                    for attr_name in models[model_name].template.attributes.keys():
                        if attr_name in revision:
                            if attr_name in existing:
                                exisiting_diffable[attr_name] = existing[attr_name]
                            revision_diffable[attr_name] = revision[attr_name]

                    self.log.info('preparing diff')

                    payload: str = "\n".join(difflib.Differ().compare(
                        json.dumps(exisiting_diffable, indent=2),
                        json.dumps(revision_diffable, indent=2),
                    ))

                    if not payload:
                        continue

                    payload_bytes = payload.encode('utf8')
                    file = io.BytesIO(payload_bytes)

                    self.log.info('starting upload')
                    for upload in metis.upload_file(
                            get_project_name(dag),
                            self.debug_bucket,
                            f"link/{dag.dag_id}/{ti.task_id}/{ti.run_id}/{model_name}/{identifier}.json",
                            file,
                            size=len(payload_bytes)
                        ):
                        self.log.info('uploading...')

class MetisEtlHelpers:
    hook: EtnaHook

    def __init__(self, tail_folders: XComArg, tail_files: XComArg, hook: EtnaHook):
        self.hook = hook
        self.tail_files = tail_files
        self.tail_folders = tail_folders

    def prepare_task_token(self, read_only: bool) -> XComArg:
        """
        Returns an XComArg that can be used to pass a task token into other tasks.
        Generally, this is used with `run_on_docker` to pass a task token into an input
        file.
        """
        type = "read" if read_only else "write"

        @task(task_id=f"prepare_{type}_task_token")
        def prepare_task_token(read_only):
            token = self.hook.get_task_auth(read_only=read_only).token.decode("utf8")
            return token

        return prepare_task_token(read_only)

    def process_and_link_matching_file(
        self,
        listed_record_matches: XComArg,
        regex: re.Pattern,
        processor: Callable[[], UpdateRequest],
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

        @link(dry_run=dry_run, task_id=processor.__name__)
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

    def link_matching_file(
        self,
        # an XComArg that lists MatchedRecordFolder and Files, usually the result of helpers.list_match_folders
        listed_record_matches: XComArg,
        # the attribute of each MatchedRecordFolder to link any matching file
        attr_name: str,
        # A regex to match against any file name inside of a MatchedRecordFolder, resulting files will be linked.
        regex: re.Pattern,
        # By default, the linking process validates with magma but does not commit.  set dry_run=False to commit.
        dry_run=True,
    ) -> XComArg:
        """
        Creates a task that iterates the provided listed_record_matches (the result of a helpers.list_match_folders
        call), searching for files whose file_name matches the given regex, linking the file path to the given
        match's record in magma at the given attr_name attribute.

        Notably, when multiple files are found, the -last- file detected will be used in magma.  If you intend to link
        several matching files, use link_matching_files for attributes that are file_collection.
        """

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
        # This regex is matched first against the folder subpath within the matched record folder before linking all matching files (including none)
        folder_path_regex: re.Pattern = re.compile(r"^$"),
        dry_run=True,
    ) -> XComArg:
        """
        Like link_matching_file, except for two key points:
        1.  all matches found are linked together as a file collection
        2.  An additional argument, folder_path_regex, is matched against the folder before linking.
            Because listed_record_matches ONLY sees a single directory's contents at a time, this
            folder_path_regex ensures the correct folder is selected before linking, in case multiple
            child directories exist for a matched record folder.
        """

        @link(attr_name, dry_run=dry_run, task_id=f"link_{attr_name}")
        def do_link(listed_matches: List[Tuple[MatchedRecordFolder, List[File]]]):
            for match, files in listed_matches:
                if folder_path_regex.match(match.match_folder_subpath):
                    matched_files = [f for f in files if file_regex.match(f.file_name)]
                    yield match, matched_files

        return do_link(listed_record_matches)

    def find_record_folders(
        self,
        # The model for which these matches represent identifiers for.
        model_name: str,
        # A regex used to match the folder path of any file or folder matched.  The part of the path that matches
        # this regex is used as the MatchedRecordFolder's root, the record name being determined by taking the last
        # path segment of the match.
        regex: re.Pattern,
        # A function that should map the matched directory name to an expected record name in magma
        # Ideally, matched folders are correct magma record names, but if some transformation is necessary,
        # this function can be specified.
        corrected_record_name: Callable[[str], str] = lambda x: x,
    ) -> XComArg:
        """
        Creates a task that produces a list of MatchedRecordFolders by finding any file or folder change that is
        found inside a folder given by the provided regex.

        This function should be called to identify folders that have a 1:1 mapping to a magma model + record combination.
        Doing so allows many linking convenience functions that can determine the record name by the containing folder's
        name.
        """

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
        """
        Creates a task that first validates the model / record name of the provided record folder matches already exist
        in magma / timur, ensuring that any incorrectly named folder does not accidentally become a magma record.
        """

        @task
        def filter_by_timur(matches):
            with self.hook.magma() as magma:
                return pickled(filter_by_exists_in_timur(magma, matches))

        return filter_by_timur(matches)

    def list_match_folders(self, matches: XComArg):
        """
        Creates a task that lists the contents of the given MatchedRecordFolders.  For instance, if a single file
        is added to a record folder, the entire contents of that record folder are listed for use by downstream
        linkers.  this is important when linking multiple files in a collection, or processing that involves
        reading in multiple files that exist in the same record folder.

        The resulting XComArg will provide a value to any consuming task a value of List[Tuple[MatchedRecordFolder, List[File]]]
        """

        @task
        def list_match_folders(matches):
            if not matches:
                return []

            project_name = matches[0].project_name
            bucket_name = matches[0].bucket_name
            with self.hook.metis() as metis:
                files_by_parent_id = { k: list(v) for k, v in itertools.groupby(metis.tail(project_name, bucket_name, 'files', folder_id=[m.folder_id for m in matches])[0], key=lambda file: file.folder_id) }
                return [
                    (m, files_by_parent_id[m.folder_id])
                    for m in matches
                ]

        return list_match_folders(matches)


def filter_by_exists_in_timur(
    magma: Magma,
    matched: List[MatchedRecordFolder],
) -> List[MatchedRecordFolder]:
    """
    Underlying function that implements helpers.filter_by_timur, takes a magma client and a list of matched record folders,
    and returns the list of matches that have a backing record in magma.
    """
    if not matched:
        return []

    project_name = matched[0].project_name
    matched_by_model_name = {}
    for m in matched:
        matched_by_model_name.setdefault(m.model_name, []).append(m)

    result: List[MatchedRecordFolder] = []

    for model_name, matches in matched_by_model_name.items():
        response = magma.retrieve(
            project_name,
            model_name=model_name,
            attribute_names="identifier",
            record_names=[m.record_name for m in matches],
        )
        result.extend(
            m for m in matches if m.record_name in response.models[model_name].documents
        )

    return result


def filter_by_record_directory(
    files_or_folders: List[Union[File, Folder]],
    directory_regex: re.Pattern,
    model_name: str,
    corrected_record_name: Callable[[str], str] = lambda x: x,
) -> List[MatchedRecordFolder]:
    """
    Underlying implementation that backs helpers.find_record_folders, searching a list of File and Folder objects
    for any whose file_path matches the directory regex, and creating MatchedRecordFolders for each unique folder_path
    found in this way.
    """
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

    tail_type: Union[Literal['files'], Literal['folders']] = 'files'
    if type == 'folder':
        tail_type = 'folders'

    files, folders = metis.tail(project_name, bucket_name, tail_type, start, end)

    if type == "file":
        return files
    else:
        return folders
