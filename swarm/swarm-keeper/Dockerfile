FROM etna-base-dev

FROM alpine
RUN mkdir /app
RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
RUN pip install -U requests[socks] yq
COPY lib/* /usr/lib/
COPY --from=0 /bin/post-to-slack.sh /usr/local/bin/post-to-slack

ENV DOCKER_VERSION 20.10.6
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/

WORKDIR /app