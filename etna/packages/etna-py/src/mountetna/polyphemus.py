import dataclasses
import typing
from io import StringIO
from functools import partial
from inspect import isgenerator
from typing import Dict, Optional, List
from serde import serialize, deserialize
from serde.json import from_json, to_json

import pandas as pd

from .etna_base import EtnaClientBase
from requests import HTTPError
from .utils.iterables import batch_iterable
from .metis import Upload, File


class Polyphemus(EtnaClientBase):
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

    def query(self, query: QueryRequest):
        response = self.session.post(
            self.prepare_url("query"),
            data=to_json(query),
            headers={"Content-Type": "application/json"},
        )
        if 'format' in query and 'tsv' == query['format']:
            pd_wrapper = partial(pd.read_csv, sep="\t")
            tsv_data = StringIO(response.text)
            return pd_wrapper(tsv_data)

        return from_json(QueryResponse, response.content)

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
