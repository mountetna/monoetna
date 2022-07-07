FROM archimedes-r-base
# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY . /app/

# RUN cd /app/archimedes-r && npm install && npm run build
RUN jetpack init && cd /app/archimedes-r

ARG FULL_BUILD=1
RUN /entrypoints/build.sh