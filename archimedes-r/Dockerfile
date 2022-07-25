# FROM development-certs
FROM etnaagent/archimedes-r-base:main

# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY archimedes-r/functions /app/archimedes-r/functions
RUN R -e "jetpack::add('dataflow', remotes='local::archimedes-r/functions/dataflow')"

# NOTE: Anything after this will almoast always get run every time
COPY . /app/