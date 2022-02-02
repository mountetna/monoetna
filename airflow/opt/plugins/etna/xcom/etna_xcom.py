from typing import Any

import pickle
from airflow.models.xcom import BaseXCom


# "Just like" normal XCom, but allows for values that execute custom logic on read, allowing for
# external references or custom logic.
class EtnaDeferredXCom:
    def execute(self):
        return None

    def __str__(self):
        return f"<{self.__class__.__name__}>"

class EtnaXCom(BaseXCom):
    @staticmethod
    def serialize_value(value: Any):
        if isinstance(value, EtnaDeferredXCom):
            return pickle.dumps(value)
        return BaseXCom.serialize_value(value)

    @staticmethod
    def deserialize_value(result: "EtnaXCom") -> Any:
        result = BaseXCom.deserialize_value(result)
        if isinstance(result, EtnaDeferredXCom):
            return result.execute()
        return result

    def orm_deserialize_value(self) -> Any:
        result = BaseXCom.deserialize_value(self)
        if isinstance(result, EtnaDeferredXCom):
            return str(result)
        return BaseXCom.deserialize_value(self)