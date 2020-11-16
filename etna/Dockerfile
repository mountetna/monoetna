FROM etna-base
ENV ETNA_GEM_DEVELOPMENT=1
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install


# RUN apt-get update &&\
#     apt-get install -y r-base \
#     libapreq2-dev \
#     apache2 \
#     apache2-dev \
#     git \
#     libpcre3 \
#     libpcre3-dev \
#     liblzma-dev \
#     bzip2 \
#     libbz2-dev \
#     zlib1g \
#     zlib1g-dev \
#     libicu-dev \
#     --no-install-recommends && \
#     rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y r-base \
    --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*
