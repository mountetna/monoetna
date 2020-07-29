set -e
export PGPASSWORD=password
echo "
CREATE DATABASE ${APP_NAME}_test
GRANT ALL PRIVILEGES ON DATABASE ${APP_NAME}_test TO developer;
" | psql -U developer -h ${APP_NAME}_db ${APP_NAME}
