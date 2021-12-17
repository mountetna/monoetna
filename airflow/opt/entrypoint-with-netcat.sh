#!/usr/bin/env bash
set -eo pipefail

rm -rf /myfifo
mkfifo /myfifo
nc -lkv 8085 </myfifo | nc -Uv /var/run/docker.sock >/myfifo &

exec /opt/airflow/entrypoint-with-secret-files.sh $@
