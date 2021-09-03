FROM apache/airflow:2.1.3

#RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
#RUN pip install -U requests[socks]
#COPY src/dags /opt/airflow/dags
#RUN pip3 install dask[complete]

COPY src/opt/airflow.cfg /opt/airflow/
COPY src/opt/webserver_config.py /opt/airflow/
COPY etna-base-dev/bin/post-to-slack.sh /usr/local/bin/post-to-slack