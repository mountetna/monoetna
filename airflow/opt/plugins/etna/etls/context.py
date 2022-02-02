from datetime import datetime
from typing import Tuple

from airflow.models.taskinstance import Context

batch_start_context_key = 'prev_data_interval_end_success'
batch_end_context_key = 'data_interval_end'

def get_batch_range(context: Context) -> Tuple[datetime, datetime]:
    start_date = context[batch_start_context_key]
    if start_date is None:
        start_date = datetime(2000, 1, 1)
    return start_date, context[batch_end_context_key]
