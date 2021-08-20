#!/usr/bin/env bash

#for replica_ip in $(getent hosts "${SENTINEL_NAME}" | awk '{print $1}')
#do
#  # skip this replica
#  test "${replica_ip}" == "${CONTAINER_IP}" && continue
#
#  # reset other sentinel
#  redis-cli -h "${replica_ip}" -p "${SENTINEL_PORT}" sentinel reset \*
#done
#

exec $@
