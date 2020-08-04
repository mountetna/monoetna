FROM library/httpd:2.4.43
WORKDIR /app

RUN apt-get update -y
RUN apt-get install -y libapache2-mod-xsendfile
COPY certs /root/certs
