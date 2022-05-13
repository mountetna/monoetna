FROM archimedes-node-base
# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY . /app/

RUN cd /app/archimedes-node && npm install && npm run build

ARG FULL_BUILD=1
RUN /entrypoints/build.sh
