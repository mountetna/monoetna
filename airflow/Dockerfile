FROM apache/airflow:2.2.2

# This drops an sqlite database in the base directory which is not used in production,
# but enables much faster testing since the migrations do not need to be run.
RUN airflow db reset --yes

# RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
USER root
RUN apt-get update
RUN apt-get install -y bash curl jq git openssh-client

USER airflow

COPY --chown=airflow:0 ./opt/requirements.txt /opt/airflow/
RUN pip3 install -r /opt/airflow/requirements.txt

COPY --chown=airflow:0 ./opt/ /opt/airflow/

ENTRYPOINT ["/opt/airflow/entrypoint-with-secret-files.sh"]
