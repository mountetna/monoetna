from typing import Mapping, Any

from airflow.utils import yaml

def etl_args(data: Mapping[str, Any], prefix="ETL__") -> Mapping[str, str]:
    env = {}

    for key, value in data.items():
        env[prefix + key.upper()] = yaml.dump(value)

    return env