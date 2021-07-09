This container starts the sample Airflow application using Celery, per the documentation here:

https://airflow.apache.org/docs/apache-airflow/stable/start/docker.html

Note that our Makefiles will also bring up and down the Airflow container.

Couple of steps to note:

1. You will want to make the `logs` container and give it the right permissions.

```
$ cd airflow
$ mkdir ./logs
$ echo -e "AIRFLOW_UID=$(id -u)\nAIRFLOW_GID=0" > .env
```

2. In the Airflow directory, run migrations and set up the database.

```
$ docker-compose up airflow-init
```

3. Now you can go back to the root `monoetna` directory and do anything like `$ make up` or `$ make -C airflow up`.
4. Edit your `/etc/hosts/` file to include `127.0.0.1 airflow.development.local`.
5. Open up a browser to `https://airflow.development.local`, and you can log in with the default user / password of `airflow` / `airflow`.
