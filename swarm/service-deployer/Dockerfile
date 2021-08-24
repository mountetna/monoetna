FROM alpine
RUN mkdir /config
RUN apk add --no-cache bash curl jq openssh-client
RUN curl -fsLo /bin/tini https://github.com/krallin/tini/releases/download/v0.19.0/tini-amd64 && \
  chmod +x /bin/tini
ENV DOCKER_VERSION 20.10.6
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/
COPY src/bin/* /usr/local/bin/
COPY etna-base-dev/bin/post-to-slack.sh /usr/local/bin/post-to-slack

WORKDIR /

CMD ["watch-and-update-service-images"]
