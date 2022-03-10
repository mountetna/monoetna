from typing import Optional

from airflow.exceptions import AirflowNotFoundException
from airflow.hooks.base import BaseHook


def find_first_valid_connection(*options: Optional[str]) -> str:
    last_e = None
    for option in options:
        if option is None:
            continue

        try:
            BaseHook.get_connection(option)
            return option
        except AirflowNotFoundException as e:
            last_e = e
    if last_e:
        raise last_e

    raise AirflowNotFoundException("No connection option present!")
