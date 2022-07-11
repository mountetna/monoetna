FROM development-certs
FROM satijalab/seurat:4.1.0

RUN apt-get update \
   && apt-get install -y \
      curl \
      libopenblas-base \
      libcurl4-openssl-dev \
      libssl-dev \
      wget \
      gnupg \
      bash-completion \
      git \
      build-essential \
      libgeos-dev \
      libmagick++-dev \
      libcairo2-dev \
      libfontconfig1-dev \
      libhdf5-dev \
      --no-install-recommends

# Add jetpack package manager and enable command line controls
RUN R -e 'install.packages("jetpack"); jetpack::cli()'

# Copy entrypoint scripts
COPY entrypoints /entrypoints
