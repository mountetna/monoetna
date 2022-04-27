FROM archimedes-node-base
# Perform these steps first to allow better caching behavior
WORKDIR /app
COPY . /app/
RUN chown -R 1000:1000 /app/node_modules
RUN npm install
ARG FULL_BUILD=1
RUN /entrypoints/build.sh
