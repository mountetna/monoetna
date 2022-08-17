FROM r-base:4.2.1

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
      libharfbuzz-dev \
      libfribidi-dev \
      --no-install-recommends

# Add jetpack package manager and enable command line controls
WORKDIR /app
COPY DESCRIPTION .
COPY renv.lock .
COPY .Rprofile .
RUN R -e 'install.packages("jetpack"); jetpack::cli()'
RUN R -e "jetpack::install()"
