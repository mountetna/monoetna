FROM etna-base
# Perform these steps first to allow better caching behavior
COPY . /app/
ARG APP_NAME
ARG FULL_BUILD=1
RUN /entrypoints/build.sh
RUN chmod +x /app/bin/magma-puma.sh && ln -sf /app/bin/magma-puma.sh /app/bin/puma.sh
