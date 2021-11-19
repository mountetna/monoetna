FROM curlimages/curl AS curl
USER root
# Tini is a little init script that handles signal processing and zombie processes, but
# otherwise should not effect system behavior.  Good for long lived docker processes.
RUN curl -fsLo /tini https://github.com/krallin/tini/releases/download/v0.19.0/tini-amd64 && \
  chmod +x /tini

FROM redis:6
COPY --from=curl /tini /usr/local/bin/
COPY bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*

ENTRYPOINT ["tini", "init-config-and-run-ha.sh"]
CMD ["redis-server"]