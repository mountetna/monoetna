from datetime import datetime, timedelta
from typing import Tuple

from airflow.models.taskinstance import Context

batch_start_context_key = 'prev_data_interval_end_success'
batch_end_context_key = 'data_interval_end'

def get_batch_range(context: Context, inclusive=True) -> Tuple[datetime, datetime]:
    start_date = context[batch_start_context_key] or datetime(2010, 1, 1)
    start_date = start_date.replace(microsecond=0)

    end_date = context[batch_end_context_key]

    # For etls where time boundaries should overlap slightly to ensure precision isn't a factor in missed
    # entries.
    if inclusive:
        end_date += timedelta(seconds=1)

    end_date = end_date.replace(microsecond=0)

    return start_date, end_date
