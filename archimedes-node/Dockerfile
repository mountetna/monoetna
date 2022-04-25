FROM archimedes-node-base
# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY . /app/
RUN npm install
ARG FULL_BUILD=1
RUN /entrypoints/build.sh
