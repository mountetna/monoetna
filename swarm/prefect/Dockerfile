FROM etna-base-dev

FROM python:3.8.3-slim
RUN pip install -U prefect prefect[redis,viz,templates]
COPY --from=0 /bin/post-to-slack.sh /usr/local/bin/post-to-slack

WORKDIR /app