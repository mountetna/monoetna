# FROM development-certs
FROM archimedes-r-base

# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY archimedes-r/functions /app/archimedes-r/functions
RUN R -e "jetpack::add('dataflow', remotes='local::archimedes-r/functions/dataflow')"
RUN R -e "jetpack::add('dataSupport', remotes='local::archimedes-r/functions/dataSupport')"

# NOTE: Anything after this will almoast always get run every time
COPY . /app/