FROM known_hosts
FROM bash_mocker

FROM apache/airflow:2.2.2-python3.8

# This drops an sqlite database in the base directory which is not used in production,
# but enables much faster testing since the migrations do not need to be run.
RUN airflow db reset --yes

# RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
USER root
# Authorize SSH Hosts here with ssh-keyscan entries
COPY --from=0 /known_hosts /root/.ssh/known_hosts
RUN chmod 0700 /root/.ssh
RUN apt-get update
RUN apt-get install -y bash curl jq git openssh-client ncat
RUN groupadd --gid 999 docker \
    && usermod -aG docker airflow

USER airflow

COPY --chown=airflow:0 --from=0 /known_hosts /opt/airflow/.ssh/known_hosts
RUN chmod 0700 /opt/airflow/.ssh
COPY --chown=airflow:0 ./opt/requirements.txt /opt/airflow/
RUN pip3 install -r /opt/airflow/requirements.txt

COPY --chown=airflow:0 ./opt/ /opt/airflow/
COPY --from=1 /mocker /usr/lib/mocker
RUN pip uninstall -y apache-airflow-providers-google apache-airflow-providers-grpc apache-airflow-providers-hashicorp apache-airflow-providers-microsoft-azure apache-airflow-providers-sendgrid apache-airflow-providers-elasticsearch apache-airflow-providers-cncf-kubernetes apache-airflow-providers-celery
RUN pip install /opt/airflow/providers/etna

ENTRYPOINT ["/opt/airflow/entrypoint-with-secret-files.sh"]
