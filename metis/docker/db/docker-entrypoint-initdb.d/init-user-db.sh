#!/bin/bash
set -e


# Add any other roles / databases here
psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
  CREATE DATABASE metis_test;
  GRANT ALL PRIVILEGES ON DATABASE metis_test TO $POSTGRES_USER;
  GRANT ALL PRIVILEGES ON DATABASE metis_development TO $POSTGRES_USER;
EOSQL
