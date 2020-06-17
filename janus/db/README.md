### Pre migration setup.

In postgres you want to run the following commands to setup the database for the migration. As you can see, we are disabling the `public` Postgres Schema and creating a new one called `private`. Make sure (if you are using Postgres) that you set the `search_path` in your `config.yml` file to the correct Postgres Schema.

```
$ sudo -i -u postgres
$ psql
postgres=# drop database janus;
postgres=# create database janus;
postgres=# \c janus;
janus=# REVOKE ALL ON schema public FROM public;
janus=# REVOKE ALL ON DATABASE janus FROM public;
janus=# GRANT CONNECT ON DATABASE janus TO developer;
janus=# CREATE SCHEMA private;
janus=# GRANT CREATE, USAGE ON SCHEMA private TO developer;
janus=# GRANT CREATE, USAGE ON SCHEMA public TO developer;
janus=# \q
$ exit
```

### Ruby Sequel migration command.

`sequel -m ./migrations postgres://[USER]:[PASS]@[HOST]/[DB NAME]?search_path=[SCHEMA]`

ex:
`sequel -m ./migrations postgres://developer:p455w0rd@localhost/janus?search_path=private`

### Post migration setup.

On the command line you want to seed your Janus DB with initial data. There is an example file here to help you get started.

WARNING! The example file contains example secret keys, you should rotate those keys first thing. If you do not know how to create new secure random strings/keys then you should not be doing this procedure. Using the keys as is from the example file is like leaving the door keys to your house in your front lock. CHANGE/ROTATE YOUR KEYS!

`$ psql -U developer -d janus < seed.sql`

### About the 'administation' project.

Administration users are necessary to the Mount Etna system. Make sure you have at least one administrator setup in the system to be able to set up the rest of the services.
