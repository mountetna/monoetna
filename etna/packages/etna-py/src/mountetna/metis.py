import dateutil

@serialize
@deserialize
@dataclasses.dataclass
class MagmaFileEntry:
    path: str = ""
    original_filename: str = ""

@serialize
@deserialize
@dataclasses.dataclass
class TailNode:
    id: int = 0
    type: str = "file"
    node_name: str = ""
    updated_at: str = ""
    parent_id: Optional[int] = None
    file_hash: Optional[str] = None
    archive_id: Optional[int] = None

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.parse(self.updated_at)

class TailResultContainer:
    bucket_name: str
    project_name: str
    nodes: typing.List[TailNode]
    parents: typing.Dict[int, TailNode]
    path_cache: typing.Dict[int, str]

    def __init__(self, bucket_name: str, project_name: str):
        self.nodes = []
        self.parents = {}
        self.bucket_name = bucket_name
        self.project_name = project_name
        self.path_cache = {}

    def add(self, node: TailNode):
        if node.type == 'parent':
            self.parents[node.id] = node
        else:
            self.nodes.append(node)

    def resolve_files(self) -> List["File"]:
        result = []
        for node in self.nodes:
            if node.type != 'file':
                continue
            result.append(File(
                file_name=node.node_name,
                file_hash=node.file_hash,
                updated_at=node.updated_at,
                file_path=self.resolve_path(node),
                folder_id=node.parent_id,
                project_name=self.project_name,
                bucket_name=self.bucket_name,
            ))

        return result

    def resolve_folders(self) -> List["Folder"]:
        result = []
        for node in self.nodes:
            if node.type != 'folder':
                continue
            result.append(Folder(
                folder_name=node.node_name,
                updated_at=node.updated_at,
                folder_path=self.resolve_path(node),
                id=node.id,
                project_name=self.project_name,
                bucket_name=self.bucket_name,
                folder_id=node.parent_id,
            ))

        return result

    def resolve_path(self, node: TailNode):
        if node.parent_id is None:
            return node.node_name

        parent_path = self.path_cache.get(node.parent_id)

        if parent_path is None:
            parent = self.parents[node.parent_id]
            parent_path = self.resolve_path(parent)
            self.path_cache[node.parent_id] = parent_path

        return f"{parent_path}/{node.node_name}"

@serialize
@deserialize
@dataclasses.dataclass
class File:
    file_name: str = ""
    file_hash: str = ""
    updated_at: str = ""
    file_path: Optional[str] = ""
    folder_id: Optional[int] = None
    project_name: str = ""
    bucket_name: str = ""
    download_url: Optional[str] = ""

    def __str__(self):
        return self.as_metis_url

    @property
    def as_magma_file_attribute(self) -> MagmaFileEntry:
        return MagmaFileEntry(
            path=self.as_metis_url,
            original_filename=self.file_name or os.path.basename(self.file_path),
        )

    @property
    def as_metis_url(self):
        return f"metis://{self.project_name}/{self.bucket_name}/{self.file_path}"

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.parse(self.updated_at)


@serialize
@deserialize
@dataclasses.dataclass
class Folder:
    id: int = 0
    folder_path: str = ""
    project_name: str = ""
    bucket_name: str = ""
    folder_name: str = ""
    updated_at: str = ""
    folder_id: int = 0

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.parse(self.updated_at)


@serialize
@deserialize
@dataclasses.dataclass
class FoldersAndFilesResponse:
    folders: List[Folder] = dataclasses.field(default_factory=list)
    files: List[File] = dataclasses.field(default_factory=list)

    def empty(self):
        return not self.folders and not self.files

    def extend(self, other: "FoldersAndFilesResponse"):
        self.folders.extend(other.folders)
        self.files.extend(other.files)


@serialize
@deserialize
@dataclasses.dataclass
class FoldersResponse:
    folders: List[Folder]


class Metis(EtnaClientBase):
    def touch_folder(
        self, project_name: str, bucket_name: str, folder_path: str
    ) -> FoldersResponse:
        response = self.session.get(
            self.prepare_url(project_name, "folder", "touch", bucket_name, folder_path)
        )
        return from_json(FoldersResponse, response.content)

    def list_folder(
        self, project_name: str, bucket_name: str, folder_path: Optional[str] = None
    ) -> FoldersAndFilesResponse:
        if folder_path:
            response = self.session.get(
                self.prepare_url(project_name, "list", bucket_name, folder_path)
            )
        else:
            response = self.session.get(
                self.prepare_url(project_name, "list", bucket_name)
            )
        return from_json(FoldersAndFilesResponse, response.content)

    def authorize_download(self, project_name: str, bucket_name: str, file_path: str) -> str:
        response = self.session.post(self.prepare_url('authorize', 'download'), json=dict(
            project_name=project_name,
            bucket_name=bucket_name,
            file_path=file_path,
        ))

        return response.json()['download_url']

    def open_file(self, file: File, binary_mode=False) -> io.BufferedReader:
        """
        See open_url, this method opens the given file's download_url using that method.
        If the given file does not have a download_url set, it will fetch one using `authorize_download`.
        """
        download_url = file.download_url
        if not download_url:
            download_url = self.authorize_download(file.project_name, file.bucket_name, file.file_path)

        return self.open_url(download_url, binary_mode)

    def open_url(self, url: str, binary_mode=False) -> io.BufferedReader:
        """
        Opens the given url for download into a context as a python io (file like) object.
        By default, when binary_mode is False, the underlying io object yields decoded strings using utf-8.
        Otherwise, the provided io object yields bytes objects.

        Note:  Ideally, this method is used in combination with 'with' syntax in python, so that the underlying
        data stream is closed after usage.  This is especially performant when code only needs to access a small
        subset of the data.

        eg:
        ```
        with metis.open_url(file.download_url) as open_file:
          for line in csv.reader(open_file):
             break
        ```
        """
        response = self.session.get(url, stream=True)
        io_obj = iterable_to_stream(response.iter_content(io.DEFAULT_BUFFER_SIZE))
        if not binary_mode:
            io_obj = io.TextIOWrapper(io_obj, encoding='utf8')
        return io_obj

    def authorize_upload(
        self, project_name: str, bucket_name: str, file_path: str
    ) -> "UploadResponse":
        response = self.session.post(
            self.prepare_url("authorize", "upload"),
            json=dict(
                project_name=project_name, bucket_name=bucket_name, file_path=file_path
            ),
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def start_upload(self, request: "UploadStartRequest") -> "UploadResponse":
        response = self.session.post(
            self.prepare_url_unsafe(request.upload_path),
            data=to_json(request),
            headers={"Content-Type": "application/json"},
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def upload_blob(self, request: "UploadBlobRequest"):
        response = self.session.post(
            self.prepare_url_unsafe(request.upload_path),
            files=encode_as_multipart(request),
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def upload_file(
        self,
        project_name: str,
        bucket_name: str,
        dest_path: str,
        file: typing.IO,
        size: Optional[int] = None,
        max_retries=3,
    ) -> typing.Iterable["Upload"]:
        if size is None:
            if file.seekable():
                file.seek(0, 2)
                size = file.tell()
                file.seek(0)
            else:
                raise ValueError(
                    "size was not specified, but file is streaming.  You must specify a size hint when using a stream."
                )

        if len(dest_path) > 1 and not self.folder_exists(project_name, bucket_name, os.path.dirname(dest_path)):
            self.create_folder(project_name, bucket_name, os.path.dirname(dest_path))
        authorization = self.authorize_upload(project_name, bucket_name, dest_path)
        upload = Upload(file=file, file_size=size, upload_path=authorization.upload_path)
        return self.upload_parts(upload, uuid.uuid4().hex, max_retries)

    def upload_parts(
        self, upload: "Upload", metis_uid: str, remaining_attempts: int, reset=False
    ) -> typing.Iterable["Upload"]:
        unsent_zero_byte_file = upload.cur == 0
        upload.resume_from(
            self.start_upload(
                UploadStartRequest(
                    upload_path=upload.upload_path,
                    file_size=upload.file_size,
                    next_blob_size=upload.next_blob_size,
                    next_blob_hash=upload.next_blob_hash,
                    metis_uid=metis_uid,
                    reset=reset,
                )
            )
        )

        while unsent_zero_byte_file or not upload.is_complete:
            try:
                current_byte_position = upload.cur
                blob_data = upload.read_next_blob()

                upload.advance_position()

                self.upload_blob(
                    UploadBlobRequest(
                        upload_path=upload.upload_path,
                        next_blob_size=upload.next_blob_size,
                        next_blob_hash=upload.next_blob_hash,
                        blob_data=blob_data,
                        metis_uid=metis_uid,
                        current_byte_position=current_byte_position,
                    )
                )

                unsent_zero_byte_file = False
            except HTTPError as e:
                if remaining_attempts > 1:
                    if e.response.status_code == 422:
                        yield from self.upload_parts(
                            upload, metis_uid, remaining_attempts - 1, True
                        )
                        return
                    else:
                        yield from self.upload_parts(
                            upload, metis_uid, remaining_attempts - 1
                        )
                        return
                raise e

            yield upload

    def create_folder(
        self, project_name: str, bucket_name: str, folder_path: str
    ) -> FoldersAndFilesResponse:
        response = self.session.post(
            self.prepare_url(project_name, "folder", "create", bucket_name, folder_path)
        )
        response_obj = from_json(FoldersAndFilesResponse, response.content)
        return response_obj

    def folder_exists(
        self, project_name: str, bucket_name: str, folder_path: str
    ) -> bool:
        try:
            self.session.get(
                self.prepare_url(project_name, "list", bucket_name, folder_path)
            )
            return True
        except HTTPError as e:
            if e.response.status_code == 422:
                return False
            raise e

    def tail(
            self,
            project_name: str,
            bucket_name: str,
            type: typing.Union[typing.Literal['folders'], typing.Literal['files']],
            batch_start: Optional[datetime] = None,
            batch_end: Optional[datetime] = None,
            folder_id: Optional[typing.Union[int, typing.List[int]]] = None,
     ) -> typing.Tuple[List[File], List[Folder]]:
        container = TailResultContainer(bucket_name, project_name)
        if batch_start and batch_end:
            args = dict(
                batch_start=batch_start.isoformat(timespec="seconds"),
                batch_end=batch_end.isoformat(timespec="seconds"),
                type=type,
            )
        else:
            args = dict(
                folder_id=folder_id,
                type=type,
            )

        response = self.session.post(self.prepare_url(project_name, 'tail', bucket_name), json=args, stream=True)
        for line in response.iter_lines():
            if line:
                container.add(from_json(TailNode, line))

        return container.resolve_files(), container.resolve_folders()

    def find(
        self,
        project_name: str,
        bucket_name: str,
        params: List[typing.Mapping],
        limit: Optional[int] = None,
        offset: Optional[int] = None,
        hide_paths=False,
    ) -> FoldersAndFilesResponse:
        limit_set = True

        if not limit:
            limit_set = False
            limit = 1000
        if not offset:
            offset = 0

        args = dict(
            params=params,
            hide_paths=hide_paths,
            limit=limit,
            offset=offset,
        )

        result = FoldersAndFilesResponse(folders=[], files=[])
        while True:
            response = self.session.post(
                self.prepare_url(project_name, "find", bucket_name), json=args
            )
            response_obj = from_json(FoldersAndFilesResponse, response.content)
            result.extend(response_obj)

            if response_obj.empty() or limit_set:
                break

            args["offset"] += max(len(response_obj.files), len(response_obj.folders))

        return result
