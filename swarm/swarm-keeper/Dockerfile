FROM alpine
RUN mkdir /app
RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
RUN pip install -U requests[socks]
COPY src/lib/* /usr/lib/
COPY etna-base-dev/bin/post-to-slack.sh /usr/local/bin/post-to-slack

WORKDIR /app