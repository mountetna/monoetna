set -e
export PGPASSWORD=password
echo "
CREATE DATABASE ${APP_NAME}_test;
GRANT ALL PRIVILEGES ON DATABASE ${APP_NAME}_test TO developer;
" | psql -U developer ${APP_NAME}_development
