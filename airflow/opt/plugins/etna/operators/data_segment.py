import dataclasses

from airflow.models.baseoperator import BaseOperator
from airflow.models.taskmixin import TaskMixin


@dataclasses.dataclass
class DagSegment(TaskMixin):
    start: BaseOperator
    end: BaseOperator

    @property
    def roots(self):
        return self.start.roots

    @property
    def leaves(self):
        return self.end.leaves

    def set_upstream(self, other):
        return self.start.set_upstream(other)

    def set_downstream(self, other):
        return self.end.set_downstream(other)

    def update_relative(self, other: "TaskMixin", upstream=True) -> None:
        if upstream:
            return self.start.update_relative(other, upstream)
        return self.end.update_relative(other, upstream)
