FROM etna-base-dev

FROM alpine

RUN apk add --no-cache bash curl git openssh-client jq

ENV DOCKER_VERSION 19.03.12

# Add the docker client
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/
COPY --from=0 /bin/post-to-slack.sh /bin/post-to-slack
RUN chmod +x /bin/*
