import dataclasses
from typing import Dict, Optional, List
from serde import serialize, deserialize
from serde.json import from_json, to_json

from .etna_base import EtnaClientBase

@serialize
@deserialize
@dataclasses.dataclass
class RulesResponse:
    rules: Dict[ str, Dict[str, str ] ] = dataclasses.field(default_factory=dict)

class Gnomon(EtnaClientBase):
    def rules(
        self, project_names: List[str] = dataclasses.field(default_factory=list)
    ):
        response = self.session.post(
            self.prepare_url("gnomon", "rules"),
            json=dict(
                'project_names' : project_names
            ),
            verify=False
        )

        return from_json(RulesResponse, response.content)
