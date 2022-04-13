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
from serde import serialize, deserialize, from_dict, to_dict
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
    match_file: Optional[File] = None
    match_folder: Optional[Folder] = None

    # serde does not directly support recursive json structures currently, this attribute
    # simply buffers the raw dictionary for the `match_parent` attribute below.
    match_parent_raw: Optional[Any] = None

    @property
    def match_parent(self) -> Optional["MatchedRecordFolder"]:
        if self.match_parent_raw is not None:
            return from_dict(MatchedRecordFolder, self.match_parent_raw)

    @match_parent.setter
    def match_parent(self, v: Any):
        if v is not None:
            self.match_parent_raw = to_dict(v)
        else:
            self.match_parent_raw = None

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
    validate_record_update: Optional[Callable[[Any, Any, List[str]], bool]]

    def __init__(
        self,
        attribute_name: Optional[str] = None,
        dry_run=True,
        hook: Optional[EtnaHook] = None,
        task_id: Optional[str] = None,
        validate_record_update: Optional[Callable[[Any, Any, List[str]], bool]] = None
    ):
        self.task_id = task_id
        self._hook = hook
        self.dry_run = dry_run
        self.attribute_name = attribute_name
        self.log = logging.getLogger("airflow.task")
        self.validate_record_update = validate_record_update

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

            updated_names = [
                f"{model_name}/{identifier}"
                for model_name, mr in value.revisions.items()
                for identifier in mr
            ]

            update.extend(value)
            result.append(match)
            self.log.info(f"Updates to {updated_names} found in {match}")

        self.log.info("Validating update...")
        self._validate_update(magma, update, models)
        self.log.info("Executing batched magma update")
        magma.update(update)
        return result

    def _validate_update(self, magma: Magma, update: UpdateRequest, models: Dict[str, Model]):
        if not self.validate_record_update:
            return

        differ = difflib.Differ()
        for model_name, revisions in update.revisions.items():
            for revision_batch in batch_iterable(revisions.keys(), 100):
                response = magma.retrieve(get_project_name(), model_name=model_name, record_names=revision_batch)
                disconnected_response = magma.retrieve(get_project_name(), model_name=model_name, record_names=revision_batch, show_disconnected=True)
                response.extend(disconnected_response)

                for identifier in revision_batch:
                    revision = update.revisions[model_name].get(identifier, {})
                    existing = response.models[model_name].documents.get(identifier, {})
                    # Compare the subset of attributes that actually exist in this revision
                    # against what might already exist for that record.
                    # Order the attributes based on the template to ensure cleaner diff as well.
                    exisiting_diffable = dict()
                    revision_diffable = dict()
                    unequal_attrs = []

                    for attr_name in models[model_name].template.attributes.keys():
                        if attr_name in revision:
                            revision_diffable[attr_name] = a = revision[attr_name]
                            if attr_name in existing:
                                exisiting_diffable[attr_name] = b = existing[attr_name]
                                if a != b:
                                    unequal_attrs.append(attr_name)
                            else:
                                unequal_attrs.append(attr_name)

                    revision_diffable['identifier'] = identifier
                    exisiting_diffable['identifier'] = identifier

                    if not self.validate_record_update(exisiting_diffable, revision_diffable, unequal_attrs):
                        diff_str = '\n'.join(differ.compare(
                            json.dumps(exisiting_diffable, indent=2).splitlines(True),
                            json.dumps(revision_diffable, indent=2).splitlines(True),
                        ))
                        raise ValueError(f"Update on {model_name}/{identifier} failed, diffset was:\n {diff_str}")


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

    def link_record_folders_to_parent(self, embedded_matches: XComArg, dry_run=True):
        """
        When record folders are embedded in metis such that parent model record folders
        exist, it is common to want to 'link' the parent model attribute based on the folder name.

        This function takes an XComArg containing a list of MatchedRecordFolder obtained via
        `find_record_folders`, *whose source arg was not None*.  In other words, consider this example:


        Imagine a directory structure like /MYSAMPLE.T1/ScRnaSeq/MYSAMPLE.T1.BLAH.BLAH
        ```
        sample_matches = helpers.find_record_folders("sample", SAMPLE_REGEX)
        with TaskGroup("sc_rna_seq_matches"):
          sc_rna_seq_matches = helpers.find_record_folders("sc_rna_seq", SC_RNA_SEQ_REGEX, source=sample_matches)
          helpers.link_record_folders_to_parent(sc_rna_esq_matches)
        ```

        In this case, FIRST task will find "/MYSAMPLE.T1" as a "sample" record folder and it will belong
        to the `XComArg` returned as `sample_matches`.

        Next, it will search under /MYSAMPLE.T1 to find the MYSAMPLE.T1.BLAH.BLAH directory as a "sc_rna_eq" model
        directory, noting the parent record folder as MYSAMPLE.T1 and that pairing will belong to sc_rna_seq_matches.

        Lastly, the `link_record_folders_to_parent` will use the association between MYSAMPLE.T1.BLAH.BLAH and MYSAMPLE.T1
        to make the association between that sc_rna_seq record and the parent sample record.
        """
        @link(dry_run=dry_run)
        def link_record_folders_to_parent(record_matches: List[MatchedRecordFolder]):
            with self.hook.magma(project_name=get_project_name()) as magma:
                magma_models = magma.retrieve(get_project_name(), hide_templates=False).models
            for record_folder_match in record_matches:
                if record_folder_match.model_name not in magma_models:
                    raise AirflowException(f"Folder #{record_folder_match} has invalid model_name")

                if record_folder_match.match_parent is None:
                    raise AirflowException(f"Folder #{record_folder_match} does not have a match parent, did you provide a source= argument to find_record_folders?")

                parent = magma_models[record_folder_match.model_name].template.parent

                yield record_folder_match, record_folder_match.as_update({parent: record_folder_match.match_parent.record_name})

        return link_record_folders_to_parent(embedded_matches)

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

        # When not provided, all the bucket's files and folders are used as the source search directories.
        # When provided, this may by a subset of files, such as an existing set of List[MatchedRecordFolder].
        # In that case, the XComArg returned by this function will include `match_parent` set.
        source: Optional[XComArg] = None
    ) -> XComArg:
        """
        Creates a task that Searches for folders that match the given a regex creates MatchedRecordFolder
        entries with the name of the directory at the far end of the match and the given model name.  The match must
        occur at a folder boundary, so for instance, r'abc/def' does not match 'abc/defg' but does match 'abc/def/c'.
        The directory matched by the last path segment of the regex is considered a 'record name' belonging to the given
        model_name, allowing other convenience methods to process and filter with that assumption.

        `source` should be used when processing inner record folders, with the outer record matches XComArg being passed.
        In that case, you should use a `with TaskGroup('---')` surrounding this task to distinguish the inner task name
        from the outer one.
        """

        if source is None:
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
        else:
            @task
            def find_record_folders(source):
                return pickled(
                    filter_by_record_directory(
                        source,
                        regex,
                        model_name=model_name,
                        corrected_record_name=corrected_record_name,
                        )
                )

            return find_record_folders(source)

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
                    (m, files_by_parent_id.get(m.folder_id) or [])
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
    files_or_folders: List[Union[File, Folder, MatchedRecordFolder]],
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
        file: Optional[File] = None
        folder: Optional[Folder] = None
        match_parent: Optional[MatchedRecordFolder] = None

        if isinstance(file_or_folder, File):
            path = os.path.dirname(file_or_folder.file_path)
            file = file_or_folder
        if isinstance(file_or_folder, Folder):
            path = file_or_folder.folder_path
            folder = file_or_folder
        if isinstance(file_or_folder, MatchedRecordFolder):
            path = file_or_folder.match_folder_subpath
            file = file_or_folder.match_file
            folder = file_or_folder.match_folder
            match_parent = file_or_folder

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
        if match_parent is not None:
            root_path = f"{match_parent.root_path}/{root_path}"

        match = MatchedRecordFolder(
            root_path=root_path,
            record_name=corrected_record_name(os.path.basename(root_path)),
            match_file=file,
            match_folder=folder,
            model_name=model_name,
        )

        match.match_parent = match_parent

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

# [2022-03-03, 14:01:37 PST] {metis.py:203} INFO - Updates to ['rna_seq/Control_Jurkat.Plate29', 'rna_seq/IPICRC041.T2.rna.live', 'rna_seq/IPICRC041.T2.rna.myeloid', 'rna_seq/IPICRC041.T2.rna.stroma', 'rna_seq/IPICRC041.T2.rna.tcell', 'rna_seq/IPICRC041.T2.rna.treg', 'rna_seq/IPICRC041.T2.rna.tumor', 'rna_seq/IPICRC041.T3.rna.live', 'rna_seq/IPICRC041.T3.rna.myeloid', 'rna_seq/IPICRC041.T3.rna.tcell', 'rna_seq/IPICRC041.T3.rna.treg', 'rna_seq/IPICRC050.T1.rna.live', 'rna_seq/IPICRC050.T1.rna.tcell', 'rna_seq/IPICRC050.T1.rna.tumor', 'rna_seq/IPICRC058.N1.rna.live', 'rna_seq/IPICRC060.N1.rna.live', 'rna_seq/IPICRC083.N1.rna.epcam', 'rna_seq/IPICRC083.N1.rna.live', 'rna_seq/IPICRC083.N1.rna.tcell', 'rna_seq/IPICRC088.N1.rna.live', 'rna_seq/IPICRC088.N1.rna.stroma', 'rna_seq/IPICRC088.N1.rna.tcell', 'rna_seq/IPICRC088.N1.rna.tumor', 'rna_seq/IPICRC094.N1.rna.live', 'rna_seq/IPICRC094.N1.rna.stroma', 'rna_seq/IPICRC094.N1.rna.tcell', 'rna_seq/IPICRC094.N1.rna.tumor', 'rna_seq/IPICRC099.N1.rna.epcam', 'rna_seq/IPICRC103.N2.rna.myeloid', 'rna_seq/IPICRC103.N2.rna.tcell', 'rna_seq/IPICRC104.T1.rna.treg', 'rna_seq/IPICRC104.T1.rna.tumor', 'rna_seq/IPICRC107.N1.rna.epcam', 'rna_seq/IPICRC107.T1.rna.treg', 'rna_seq/IPICRC109.N2.rna.live', 'rna_seq/IPICRC109.N2.rna.tcell', 'rna_seq/IPICRC109.N2.rna.tumor', 'rna_seq/IPICRC109.T2.rna.live', 'rna_seq/IPICRC109.T2.rna.tcell', 'rna_seq/IPICRC111.N1.rna.tcell', 'rna_seq/IPICRC111.T1.rna.stroma', 'rna_seq/IPICRC112.T2.rna.myeloid', 'rna_seq/IPICRC112.T2.rna.tcell', 'rna_seq/IPICRC116.N2.rna.myeloid', 'rna_seq/IPICRC116.N2.rna.tumor', 'rna_seq/IPICRC116.T2.rna.live', 'rna_seq/IPICRC117.N1.rna.live', 'rna_seq/IPICRC117.N1.rna.myeloid', 'rna_seq/IPICRC117.N1.rna.tcell', 'rna_seq/IPICRC117.N1.rna.tumor', 'rna_seq/IPICRC117.T1.rna.live', 'rna_seq/IPICRC117.T1.rna.myeloid', 'rna_seq/IPICRC117.T1.rna.tcell', 'rna_seq/IPICRC120.T2.rna.tcell', 'rna_seq/IPICRC125.T1.rna.treg', 'rna_seq/IPIHNSC066.T1.rna.live', 'rna_seq/IPIHNSC066.T1.rna.stroma', 'rna_seq/IPIHNSC066.T1.rna.tcell', 'rna_seq/IPIHNSC066.T1.rna.tumor', 'rna_seq/IPIHNSC091.T1.rna.live', 'rna_seq/IPIHNSC091.T1.rna.myeloid', 'rna_seq/IPIHNSC091.T1.rna.tcell', 'rna_seq/IPIHNSC091.T1.rna.treg', 'rna_seq/IPIHNSC093.T1.rna.myeloid', 'rna_seq/IPIHNSC093.T1.rna.treg', 'rna_seq/IPIHNSC093.T1.rna.tumor', 'rna_seq/IPIHNSC095.T1.rna.live', 'rna_seq/IPIHNSC095.T1.rna.tcell', 'rna_seq/IPIHNSC097.T1.rna.live', 'rna_seq/IPIHNSC097.T1.rna.myeloid', 'rna_seq/IPIHNSC097.T1.rna.treg', 'rna_seq/IPIHNSC102.T1.rna.live', 'rna_seq/IPIHNSC102.T1.rna.myeloid', 'rna_seq/IPIHNSC102.T1.rna.tcell', 'rna_seq/IPIHNSC102.T1.rna.tumor', 'rna_seq/IPIMEL317.T1.rna.tcell', 'rna_seq/IPIMEL317.T1.rna.tumor', 'rna_seq/IPIMEL323.N1.rna.live', 'rna_seq/IPIMEL323.N1.rna.myeloid', 'rna_seq/IPIMEL323.N1.rna.tcell', 'rna_seq/IPIMEL337.T2.rna.myeloid', 'rna_seq/IPIMELB004.T2.rna.live', 'rna_seq/IPIMELB004.T2.rna.myeloid', 'rna_seq/IPIMELB004.T2.rna.tcell', 'rna_seq/IPIMELB004.T2.rna.treg', 'rna_seq/IPIMELB004.T2.rna.tumor', 'rna_seq/IPIPDAC050.T1.rna.tcell', 'rna_seq/IPISRC028.T1.rna.stroma', 'rna_seq/IPISRC028.T1.rna.tcell', 'rna_seq/IPISRC030.T1.rna.myeloid', 'rna_seq/IPISRC030.T1.rna.stroma', 'rna_seq/IPISRC030.T1.rna.tcell', 'rna_seq/IPISRC030.T1.rna.tumor', 'rna_seq/IPISRC032.T1.rna.live'] found in metis://ipi/data/bulkRNASeq/processed/plate29_rnaseq_new/results/rnaseq_table.tsv