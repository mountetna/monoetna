FROM archimedes-r-base
# Perform these steps first to allow better caching behavior
WORKDIR /app
RUN jetpack init
COPY renv.lock /app
COPY DESCRIPTION /app
COPY .Rprofile /app
COPY archimedes-r/functions /app/archimedes-r/functions
RUN jetpack install

# NOTE: Anything after this will almoast always get run every time
COPY . /app/

ARG FULL_BUILD=1
RUN /entrypoints/build.sh