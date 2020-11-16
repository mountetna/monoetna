FROM etna-base
ENV ETNA_GEM_DEVELOPMENT=1
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install

RUN apt-get update && \
    apt-get install -y r-base \
    --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*
