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

@serialize
@deserialize
@dataclasses.dataclass
class ListEtlConfigsResponse:
    configs: List[Dict] = dataclasses.field(default_factory=list)

class Polyphemus(EtnaClientBase):
    def list_all_etl_configs(
        self, job_type: Optional[str] = None
    ):
        args = dict()

        if job_type is not None:
            args["job_type"] = job_type
        response = self.session.post(
            self.prepare_url("api", "etl", "configs"),
            json=args
        )

        return from_json(ListEtlConfigsResponse, response.content)
