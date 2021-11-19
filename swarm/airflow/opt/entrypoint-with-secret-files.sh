#!/usr/bin/env bash
set -eo pipefail

if [ -f /run/secrets/airflow-postgres-password ]; then
  export AIRFLOW__CORE__SQL_ALCHEMY_CONN="postgresql://developer:$(cat /run/secrets/airflow-postgres-password)@airflow_postgres:5432/airflow"
  export AIRFLOW__CORE__EXECUTOR=LocalExecutor
fi

if [ -f /run/secrets/airflow-fernet-key ]; then
  export AIRFLOW__CORE__FERNET_KEY=$(cat /run/secrets/airflow-fernet-key)
fi

echo $@
exec /entrypoint $@
