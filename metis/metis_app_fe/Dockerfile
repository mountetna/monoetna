FROM etna-apache
ARG APP_NAME
COPY . /app/build
RUN curl -L https://github.com/mikefarah/yq/releases/download/v4.53.3/yq_linux_amd64.tar.gz | tar xz && mv yq_linux_amd64 /usr/local/bin/yq
COPY entrypoints /entrypoints
RUN /entrypoints/build_app_fe.sh
