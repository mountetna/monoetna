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
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version = "0.17.3", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("cli", version = "3.6.1", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("R6", version = "2.5.1", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("rprojroot", version = "2.0.3", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("desc", version = "1.4.2", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("docopt", version = "0.7.1", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
RUN R -e 'remotes::install_version("jetpack", version = "0.5.5", repos = "http://cran.us.r-project.org", upgrade = FALSE)'
COPY DESCRIPTION .
COPY renv.lock .
COPY .Rprofile .
RUN R -e 'jetpack::cli(); jetpack::install()'