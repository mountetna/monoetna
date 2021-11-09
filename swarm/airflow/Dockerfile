FROM etna-base-dev

FROM apache/airflow:2.1.3

#RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
#RUN pip install -U requests[socks]
#COPY src/dags /opt/airflow/dags
# RUN pip3 install dask[complete]
RUN pip3 install docker pytest airflow-code-editor black

COPY --chown=airflow src/opt/ /opt/airflow/
# RUN chown -R airflow /opt/airflow/
COPY --from=0 /bin/post-to-slack.sh /usr/local/bin/post-to-slack
ENTRYPOINT ["/opt/airflow/entrypoint-with-secret-files.sh"]
