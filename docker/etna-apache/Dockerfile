FROM library/httpd:2.4.43
WORKDIR /app

ENV DOCKERIZE_VERSION 0.5.0
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz

RUN apt-get update -y
RUN apt-get install -y libapache2-mod-xsendfile libapache2-mod-shib2 curl
RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )

COPY certs /root/certs
COPY src/httpd.conf /usr/local/apache2/conf/httpd.conf
RUN ln -sf /usr/opt/shibboleth/shibboleth2.xml /etc/shibboleth/shibboleth2.xml
RUN ln -sf /usr/opt/shibboleth/incommon.pem /etc/shibboleth/incommon.pem
RUN ln -sf /usr/opt/shibboleth/sp-cert.pem /etc/shibboleth/sp-cert.pem
RUN ln -sf /usr/opt/shibboleth/sp-key.pem /etc/shibboleth/sp-key.pem
RUN ln -sf /usr/opt/shibboleth/idp-metadata.xml /etc/shibboleth/idp-metadata.xml
RUN ln -sf /usr/opt/shibboleth/attribute-map.xml /etc/shibboleth/attribute-map.xml
