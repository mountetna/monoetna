from typing import Any, TypeVar, cast, Generic

import pickle

from airflow.models.xcom import BaseXCom


# "Just like" normal XCom, but allows for values that execute custom logic on read, allowing for
# external references or custom logic.
from serde.json import to_json


class EtnaXComValue:
    def execute(self):
        return None

    def __str__(self):
        return f"<{self.__class__.__name__}>"

T = TypeVar("T")

def pickled(v: T) -> T:
    return cast(T, _Pickled(v))

# Explicit pickling
class _Pickled(EtnaXComValue):
    def __init__(self, value):
        self.value = value

    def execute(self):
        return self.value

    def __str__(self):
        if isinstance(self.value, list):
            return f"{len(self.value)} result(s)"

class EtnaXCom(BaseXCom):
    @staticmethod
    def serialize_value(value: Any):
        if isinstance(value, EtnaXComValue):
            return pickle.dumps(value)
        # Enables dataclass serialization
        return to_json(value).encode('UTF-8')

    @staticmethod
    def deserialize_value(result: "EtnaXCom") -> Any:
        # orm_deserialize_value deferred loading of pickled object.
        if isinstance(result.value, str):
            result.value = result._original_value
        # Already been deserialized.
        if not isinstance(result.value, bytes):
            return result.value

        result = BaseXCom.deserialize_value(result)
        if isinstance(result, EtnaXComValue):
            return result.execute()
        return result

    def orm_deserialize_value(self) -> Any:
        result = BaseXCom.deserialize_value(self)
        if isinstance(result, EtnaXComValue):
            self._original_value = self.value
            return str(result)
        return BaseXCom.deserialize_value(self)