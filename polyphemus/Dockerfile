FROM etna-base

# Install argo client
RUN curl -sLO https://github.com/argoproj/argo-workflows/releases/download/v3.4.8/argo-linux-amd64.gz && \
    gunzip argo-linux-amd64.gz && \
    chmod +x argo-linux-amd64 && \
    mv ./argo-linux-amd64 /usr/local/bin/argo

# Perform these steps first to allow better caching behavior
COPY . /app/
ARG APP_NAME
ARG FULL_BUILD=1
RUN /entrypoints/build.sh
RUN npm run build
