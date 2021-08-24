FROM alpine
RUN apk add --no-cache bash diffutils curl
RUN mkdir /tests
RUN curl -fsLo /tests/bach.sh https://github.com/bach-sh/bach/raw/master/bach.sh
RUN chmod +r /tests/bach.sh
RUN chmod +x /tests

WORKDIR /tests
CMD ["bash", "/tests/main.sh"]
