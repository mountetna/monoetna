FROM r-base:4.2.2

RUN apt-get update \
   && apt-get install -y \
      curl \
      libcurl4-openssl-dev \
      libssl-dev \
      libopenblas-dev \
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

WORKDIR /app
# Add jetpack package manager, enable command line controls, and instll packages based on past package uses.
COPY DESCRIPTION .
COPY renv.lock .
COPY .Rprofile .
RUN R -e 'install.packages("jetpack"); jetpack::cli()'
RUN R -e "jetpack::install()"
