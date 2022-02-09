from datetime import datetime, timedelta
from typing import Tuple, Optional

import dateutil.parser
from airflow.models import Variable
from airflow.models.taskinstance import Context, TaskInstance
from airflow.utils.timezone import utc

