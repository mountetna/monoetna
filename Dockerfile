FROM etna-base
ENV ETNA_GEM_DEVELOPMENT=1
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install
