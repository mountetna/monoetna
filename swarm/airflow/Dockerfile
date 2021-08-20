FROM curlimages/curl AS curl
USER root
RUN curl -fsLo /tini https://github.com/krallin/tini/releases/download/v0.19.0/tini-amd64 && \
  chmod +x /tini

# Change this to upgrade apache airflow.
FROM apache/airflow:2.1.2-python3.9
COPY --from=curl /tini /usr/local/bin/

COPY src/requirements.txt /opt/airflow/requirements.txt
RUN pip install -r /opt/airflow/requirements.txt

COPY src/entrypoints/init-sentinels-and-run.sh /usr/local/bin/

ENTRYPOINT ["tini", "init-sentinels-and-run.sh"]
