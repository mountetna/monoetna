FROM development-certs

FROM etna-apache
COPY development-edge-httpd.conf /usr/opt/httpd.conf.d/main.conf
COPY --from=0 /certs /root/certs
