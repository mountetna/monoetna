#!/usr/bin/env bash
# Borrowed and adjusted from https://github.com/nickadam/redis-ha/blob/main/docker-entrypoint.sh

SENTINEL_CONFIG="/data/sentinel.conf"
REDIS_CONFIG="/data/redis.conf"

PASSWORD=$(cat "${PASSWORD_FILE}")
CONTAINER_IP=$(getent hosts ${HOSTNAME} | awk '{print $1}')

# takes primary IP as argument and connects to it
create_sentinel_config(){
cat <<EOF > "${SENTINEL_CONFIG}"
sentinel monitor ${INSTANCE_NAME} ${1} "${REDIS_PORT}" "${QUORUM}"
sentinel down-after-milliseconds ${INSTANCE_NAME} 5000
sentinel failover-timeout ${INSTANCE_NAME} 60000
sentinel parallel-syncs ${INSTANCE_NAME} 1
sentinel auth-pass ${INSTANCE_NAME} ${PASSWORD}
EOF
}

create_redis_config(){
cat <<EOF > "${REDIS_CONFIG}"
dir /data
appendonly yes
masterauth ${PASSWORD}
user default on +@all ~* >${PASSWORD}
EOF
}

echo "starting up..."

# we are starting a redis server
if [ "$1" == "redis-server" ]
then
  echo "Creating redis config"
  create_redis_config

  # wait a while for the other replicas to start
  sleep 10

  while true
  do
    echo "awaiting primary or checking for low ip"
    echo "===="
    echo "Primary" $(get-primary)
    echo "Replicas" $(getent hosts redis)
    echo "===="
    # check if there is a primary running and connect to it
    primary=$(get-primary)
    test ! -z "${primary}" && \
    echo "Starting redis server ${CONTAINER_IP} as replica of ${primary}"  && \
    exec redis-server "${REDIS_CONFIG}" --replicaof "${primary}" "${REDIS_PORT}"

    # couldn't find primary start the lowest IP as the primary
    low_ip=$(getent hosts redis | awk '{print $1}' | sort -n | head -n 1)
    # this is the low ip, start as primary
    test "${low_ip}" == "${CONTAINER_IP}" && \
    echo "Starting redis server ${CONTAINER_IP}" && \
    exec redis-server "${REDIS_CONFIG}"

    # wait a little while and try again
    sleep 3
  done
fi

# we are starting a redis sentinel server
if [ "$1" == "redis-sentinel" ]
then
  # wait a while for the other replicas to start
  sleep 10

  while true
  do
    # find the primary and create config for it
    primary=$(get-primary)
    test ! -z "${primary}" && \
    echo "Starting redis sentinel ${CONTAINER_IP} monitoring ${primary}"  && \
    create_sentinel_config "${primary}" && \
    break

    # wait a little while and try again
    sleep 3
  done

  # clean up any dead sentinels by resetting others
  for replica_ip in $(getent hosts "${SENTINEL_NAME}" | awk '{print $1}')
  do
    # skip this replica
    test "${replica_ip}" == "${CONTAINER_IP}" && continue

    # reset other sentinel
    redis-cli -h "${replica_ip}" -p "${SENTINEL_PORT}" sentinel reset \*
  done

  # start redis-sentinel
  exec redis-sentinel "${SENTINEL_CONFIG}"
fi
