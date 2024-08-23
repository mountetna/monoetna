FROM bash_mocker
FROM etna-base-dev

FROM alpine
RUN mkdir /app
RUN apk add --no-cache bash curl jq py-setuptools py-pip git openssh-client vim
RUN pip install -U --break-system-packages requests[socks] yq rich httpie
COPY lib/* /usr/lib/
COPY --from=0 /mocker /usr/lib/mocker
COPY --from=1 /bin/post-to-slack.sh /usr/local/bin/post-to-slack

ENV DOCKER_VERSION 20.10.6
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/

WORKDIR /app