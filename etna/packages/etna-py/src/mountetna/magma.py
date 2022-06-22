import dataclasses
import typing
from typing import Dict, Optional, List
from serde import serialize, deserialize
from serde.json import from_json, to_json
from .etna_base import EtnaClientBase
from requests import HTTPError
from .utils.iterables import batch_iterable
from .metis import Upload

@serialize
@deserialize
@dataclasses.dataclass
class Attribute:
    attribute_type: str = ""
    link_model_name: Optional[str] = None
    match: Optional[str] = None
    restricted: Optional[bool] = None
    read_only: Optional[bool] = None
    hidden: Optional[bool] = None
    validation: typing.Any = None


@serialize
@deserialize
@dataclasses.dataclass
class Template:
    attributes: Dict[str, Attribute]
    name: str = ""
    identifier: str = ""
    version: int = 0
    parent: str = ""


@serialize
@deserialize
@dataclasses.dataclass
class Model:
    documents: Dict[str, Dict[str, typing.Any]] = dataclasses.field(default_factory=dict)
    template: Optional[Template] = None
    count: int = 0

    def empty(self):
        return len(self.documents) == 0

    def extend(self, other: "Model"):
        if other.template:
            self.template = other.template
        self.count = other.count
        self.documents.update(other.documents)


@serialize
@deserialize
@dataclasses.dataclass
class UpdateRequest:
    revisions: Dict[str, Dict[str, Dict[str, typing.Any]]] = dataclasses.field(
        default_factory=dict
    )
    project_name: str = ""
    dry_run: bool = False

    def __bool__(self):
        return not self.empty()

    def shallow_copy(self, revisions: Optional[Dict[str, Dict[str, Dict[str, typing.Any]]]] = None) -> "UpdateRequest":
        if revisions is None:
            revisions = self.revisions

        return UpdateRequest(project_name=self.project_name, dry_run=self.dry_run, revisions=dict(**revisions))

    def includes(self, model_name: str, record_name: str):
        return record_name in self.revisions.get(model_name, {})

    def __iter__(self) -> typing.Iterator[typing.Tuple[str, str]]:
        for model_name, revisions in self.revisions.items():
            for record_name in revisions.keys():
                yield model_name, record_name

    def sample_revision_tree(self, model_name: str, record_name: str, models: Dict[str, Model], acc: "UpdateRequest", destructive=False):
        """
        Attempts to sample the complete graph of related objects reference by the given model_name / record_name in this
        update, copying it and all its referenced links to acc.

        When destructive=True is set, it also removes all copied models from this update.
        """
        if not self.includes(model_name, record_name):
            return

        if model_name not in models:
            return

        if acc.includes(model_name, record_name):
            return

        # Ensure a record to prevent recursive re-entry
        acc.update_record(model_name, record_name)

        model = models[model_name]

        if destructive:
            record = self.revisions[model_name].pop(record_name)
        else:
            record = self.update_record(model_name, record_name)

        # copy table values into the acc
        for attr_name, attr in model.template.attributes.items():
            if attr_name not in record:
                continue
            val = record[attr_name]

            if attr.link_model_name and attr.attribute_type != "table":
                if isinstance(val, list):
                    for child_id in val:
                        self.sample_revision_tree(attr.link_model_name, child_id, models, acc)
                if isinstance(val, str):
                    self.sample_revision_tree(attr.link_model_name, val, models, acc)
            elif attr.link_model_name and attr.attribute_type == "table":
                if isinstance(val, list):
                    for table_id in val:
                        acc.append_table(
                            model_name,
                            record_name,
                            attr.link_model_name,
                            self.update_record(attr.link_model_name, table_id),
                            attr_name
                        )

    def __len__(self):
        return sum(len(rev) for rev in self.revisions.values())

    def empty(self):
        return len(self) == 0

    def extend(self, other: "UpdateRequest"):
        for model_name, docs in other.revisions.items():
            for record_name, revision in docs.items():
                self.update_record(model_name, record_name, revision)

    # TODO: Add more basic validations?
    def validate(self, models: Dict[str, Model]) -> typing.Iterable[str]:
        for model_name, docs in self.revisions.items():
            if model_name not in models:
                yield f"Model '{model_name}' does not exist"
                continue

            template = models[model_name].template
            for record_name, revision in docs.items():
                for attr, value in revision.items():
                    if attr not in template.attributes:
                        yield f"Model '{model_name}', Attribute '{attr}' does not exist"
                        continue

    def update_record(
        self, model_name: str, record_name: str, attrs: Dict[str, typing.Any] = dict()
    ) -> Dict[str, typing.Any]:
        attrs = {k: normalize_for_magma(v) for k, v in attrs.items()}
        record = self.revisions.setdefault(model_name, {}).setdefault(record_name, {})
        record.update(**attrs)
        return record

    def append_table(
        self,
        parent_model_name: str,
        parent_record_name: str,
        model_name: str,
        attrs: Dict[str, typing.Any],
        attribute_name: str,
    ):
        parent_revision = self.update_record(parent_model_name, parent_record_name, {})
        table = parent_revision.setdefault(attribute_name, [])
        id = f"::{model_name}{len(self.revisions.setdefault(model_name, {})) + 1}"
        table.append(id)
        self.update_record(model_name, id, attrs)
        return id

@serialize
@deserialize
@dataclasses.dataclass
class RetrievalResponse:
    models: Dict[str, Model] = dataclasses.field(default_factory=dict)

    def empty(self):
        for model in self.models.values():
            if not model.empty():
                return False
        return True

    def extend(self, other: "RetrievalResponse"):
        for model_name, model in other.models.items():
            self.models.setdefault(model_name, Model()).extend(model)

class Magma(EtnaClientBase):
    def retrieve(
        self,
        project_name: str,
        model_name="all",
        attribute_names="all",
        record_names: typing.Union[List[str], "all"] = [],
        page: Optional[int] = None,
        page_size: Optional[int] = None,
        order: Optional[str] = None,
        filter: Optional[str] = None,
        hide_templates=True,
        show_disconnected=False,
    ) -> RetrievalResponse:
        args = dict(
            project_name=project_name,
            model_name=model_name,
            attribute_names=attribute_names,
            hide_templates=hide_templates,
            record_names=record_names,
            show_disconnected=show_disconnected,
        )

        if page is None and record_names:
            page = 1
        if page_size is None and record_names:
            page_size = self.batch_size

        if page is not None:
            args["page"] = page
        if page_size is not None:
            args["page_size"] = page_size

        if order is not None:
            args["order"] = order
        if filter is not None:
            args["filter"] = filter

        response = RetrievalResponse(models={})

        def consume_pages():
            while True:
                try:
                    r = self.session.post(self.prepare_url("retrieve"), json=args)
                except HTTPError as e:
                    if b"not found" in e.response.content:
                        if e.response.status_code == 422:
                            break
                    raise e
                paged = from_json(RetrievalResponse, r.content)
                response.extend(paged)

                if "page" not in args:
                    break

                if paged.empty():
                    break
                args["page"] += 1

        if record_names and isinstance(record_names, list):
            for record_name_batch in batch_iterable(record_names, self.batch_size):
                args["record_names"] = record_name_batch
                args["page"] = 1
                consume_pages()
        else:
            args["record_names"] = record_names
            consume_pages()

        return response

    batch_size = 300

    def update(self, update: UpdateRequest):
        response = self.session.post(
            self.prepare_url("update"),
            data=to_json(update),
            headers={"Content-Type": "application/json"},
        )

        return from_json(RetrievalResponse, response.content)


def normalize_for_magma(value):
    # Exhaust generators (such as uploads) before consuming their last yielded value
    if isgenerator(value):
        for value in value:
            pass

    if isinstance(value, File):
        return value.as_magma_file_attribute

    if isinstance(value, Upload):
        magma_file_attribute = value.as_magma_file_attribute
        if magma_file_attribute is None:
            raise ValueError(
                f"Could not determine metis path for upload: '{value.upload_path}', possible bug in etna client."
            )
        return magma_file_attribute

    if isinstance(value, Upload):
        return

    if isinstance(value, list):
        return [normalize_for_magma(v) for v in value]

    return value
