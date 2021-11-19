FROM etna-base-dev

FROM alpine
RUN mkdir /app
RUN apk add --no-cache bash curl jq httpie py-setuptools py-pip git openssh-client
RUN pip install -U requests[socks]
COPY lib/* /usr/lib/
COPY --from=0 /bin/post-to-slack.sh /usr/local/bin/post-to-slack

WORKDIR /app