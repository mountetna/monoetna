#!/bin/bash
while ! mysqladmin ping --silent; do
  sleep 1
done

mysql -e "CREATE DATABASE IF NOT EXISTS slurm_acct_db;"

# Create the user if it doesn't exist and grant privileges
mysql -e "CREATE USER IF NOT EXISTS 'slurm'@'localhost' IDENTIFIED BY 'password';"
mysql -e "GRANT ALL PRIVILEGES ON slurm_acct_db.* TO 'slurm'@'localhost';"

# Apply the changes
mysql -e "FLUSH PRIVILEGES;"