"""
Improves upon the base airflow xcom with
1. lzma compression of xcom values
2. Ability to subclass EtnaXComValue to support
    a. custom string summary in UI
    b. custom deferred value expansion via execute
"""
import logging
from typing import Any, TypeVar, cast

import pickle
import lzma

from airflow.models.xcom import BaseXCom

class EtnaXComValue:
    def execute(self):
        return None

    def __str__(self):
        return f"<{self.__class__.__name__}>"


T = TypeVar("T")


def pickled(v: T) -> T:
    if isinstance(v, EtnaXComValue):
        return v
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
        return str(self.value)


class EtnaXCom(BaseXCom):
    @staticmethod
    def serialize_value(value: Any):
        if not isinstance(value, EtnaXComValue):
            value = pickled(value)
        log = logging.getLogger('airflow.task')
        log.info('Compressing pickled value...')
        compressed = lzma.compress(pickle.dumps(value))
        log.info('Compression complete.')
        return compressed

    @staticmethod
    def deserialize_value(result: "EtnaXCom") -> Any:
        # orm_deserialize_value deferred loading of pickled object.
        if isinstance(result.value, str) and hasattr(result, '_original_value'):
            result.value = result._original_value
        # Already been deserialized.
        if not isinstance(result.value, bytes):
            return result.value

        result = pickle.loads(lzma.decompress(result.value))
        if isinstance(result, EtnaXComValue):
            return result.execute()
        return result

    def orm_deserialize_value(self) -> Any:
        result = pickle.loads(lzma.decompress(self.value))
        if isinstance(result, EtnaXComValue):
            self._original_value = self.value
            return str(result)
        return self.value
