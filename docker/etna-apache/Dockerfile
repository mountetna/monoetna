FROM library/httpd:2.4.43
WORKDIR /app

ENV DOCKERIZE_VERSION 0.5.0
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz

RUN apt-get update -y
RUN apt-get install -y libapache2-mod-xsendfile libapache2-mod-shib2 curl
RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )

COPY certs /root/certs
COPY httpd.conf /usr/local/apache2/conf/httpd.conf

