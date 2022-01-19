#!/usr/bin/env bash
set -eo pipefail

if [ -f /run/secrets/airflow-postgres-password ]; then
  export AIRFLOW__CORE__SQL_ALCHEMY_CONN="postgresql://developer:$(cat /run/secrets/airflow-postgres-password)@airflow_postgres:5432/airflow"
#  export AIRFLOW__CORE__EXECUTOR=LocalExecutor
fi

if [ -f /run/secrets/airflow-fernet-key ]; then
  export AIRFLOW__CORE__FERNET_KEY=$(cat /run/secrets/airflow-fernet-key)
fi

if [ -n "$USE_MOCKS" ]; then
  echo "Preparing to use mocks"
  source /usr/lib/mocker
  mock git
#  export RECORD_CENSURE="${AIRFLOW_GIT_TEST_PASSWORD:-xxx}"
#  if [ -f "$AIRFLOW_GIT_TEST_PK_FILE" ]; then
#    export RECORD_CENSURE="${}"
#  fi
  loadRecording /opt/airflow/plugins/etna/tests/pytest.recording
  echo $PATH
fi


echo $@
exec /entrypoint $@
