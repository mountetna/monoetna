import logging
import re
from typing import List


from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context

from etna.etls.etl_task_batching import get_batch_range
from etna.hooks.cat import RsyncEntry, Cat, CatHook
from etna.etls.remote_helpers_base import RemoteHelpersBase


class CatEtlHelpers(RemoteHelpersBase):
    hook: CatHook

    def __init__(self, hook: CatHook):
        self.hook = hook
        self.log = logging.getLogger("airflow.task")


def load_cat_files_batch(
    cat: Cat,
    magic_string: re.Pattern
) -> List[RsyncEntry]:
    return _load_cat_files_batch(cat, magic_string)


def _load_cat_files_batch(
    cat: Cat,
    magic_string: re.Pattern
) -> List[RsyncEntry]:
    context: Context = get_current_context()
    _, end = get_batch_range(context)

    log = logging.getLogger("airflow.task")
    log.info(
        f"Searching for CAT data that has not been ingested as of {end.isoformat(timespec='seconds')}"
    )

    files = cat.tail(magic_string=magic_string)

    return files
