#!/usr/bin/env bash
set -eo pipefail

# Unfortunately, the python VCR library has trouble recording http over unix socket requests.
# This entrypoint, which can be substituted for the default one, simply sets up (with netcat)
# a tcp proxy to the underlying docker socket, which the tests are setup to consume.

rm -rf /myfifo
mkfifo /myfifo
nc -lkv 8085 </myfifo | nc -Uv /var/run/docker.sock >/myfifo &

exec /opt/airflow/entrypoint-with-secret-files.sh $@
