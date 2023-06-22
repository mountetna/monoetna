from dataclasses import dataclass, field
import typing
from datetime import datetime, timedelta
from io import StringIO
from functools import partial
from inspect import isgenerator
from typing import Dict, Optional, List
from serde import serialize, deserialize
from serde.json import from_json

import pandas as pd

from .etna_base import EtnaClientBase
from requests import HTTPError
from .utils.iterables import batch_iterable
from .metis import Upload, File

@serialize
@deserialize
@dataclass
class EtlConfigResponse:
    ran_at: Optional[str] = None
    status: str = ""
    created_at: str = ""
    comment: str = ""
    project_name: str = ""
    name: str = ""
    etl: str = ""
    config_id: int = 0
    version_number: int = 0
    config: Dict = field(default_factory=dict)
    secrets: Dict = field(default_factory=dict)
    params: Dict = field(default_factory=dict)
    run_interval: int = -1

    def should_run(self, start, end):
        # run_never and run_archive
        if self.run_interval < 0:
            return False

        # run_once always run
        if self.run_interval == 0:
            return True

        # new intervals always run
        if not self.ran_at:
            return True

        ran_at = datetime.fromisoformat(self.ran_at)

        if ran_at + timedelta(0, self.run_interval) < start:
            return True

        return False

@serialize
@deserialize
@dataclass
class ListEtlConfigsResponse:
    configs: List[EtlConfigResponse] = field(default_factory=list)

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

    def etl_update(
        self,
        project_name: str,
        config_id: str,
        **update
    ):
        response = self.session.post(
            self.prepare_url("api", "etl", project_name, "update", str(config_id)),
            json=update
        )

        return from_json(EtlConfigResponse, response.content)

    def etl_add_output(
        self,
        project_name: str,
        config_id: str,
        append=False,
        output: str =""
    ):
        response = self.session.post(
            self.prepare_url("api", "etl", project_name, "output", str(config_id)),
            json=dict(append=append, output=output)
        )

        return from_json(EtlConfigResponse, response.content)
