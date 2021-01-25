FROM library/httpd:2.4.43
WORKDIR /app

ENV DOCKERIZE_VERSION 0.5.0
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz

RUN apt-get update -y
RUN apt-get install -y libapache2-mod-xsendfile libapache2-mod-shib2 curl
RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )

COPY certs /root/certs
COPY src/httpd.conf /usr/local/apache2/conf/httpd.conf

# In production, shibboleth environemnt is mounted with a bind mount into the /usr/opt directory.
# We can't just mount that shibboleth prod directly into /etc/shibboleth because we only want to inject
# a subset of configuration, not replacing the entire thing.
# The solution is to have a separate directory that is mounted, but links into the individual files
# we want injected.  BUT we cannot optionally configure shibboleth in Apache due to its configuration
# language not supporting lazy evaluation, so we are ALSO forced to leave a plausible shibboleth2.xml
# in case we are not backing the link (which has to be created at a build time, which is common
# to all environments).  Yeesh.
RUN mkdir -p /usr/opt/shibboleth/
RUN cp /etc/shibboleth/shibboleth2.xml /usr/opt/shibboleth/shibboleth2.xml
RUN ln -sf /usr/opt/shibboleth/shibboleth2.xml /etc/shibboleth/shibboleth2.xml
RUN ln -sf /usr/opt/shibboleth/incommon.pem /etc/shibboleth/incommon.pem
RUN ln -sf /usr/opt/shibboleth/sp-cert.pem /etc/shibboleth/sp-cert.pem
RUN ln -sf /usr/opt/shibboleth/sp-key.pem /etc/shibboleth/sp-key.pem
RUN ln -sf /usr/opt/shibboleth/idp-metadata.xml /etc/shibboleth/idp-metadata.xml
RUN ln -sf /usr/opt/shibboleth/attribute-map.xml /etc/shibboleth/attribute-map.xml
