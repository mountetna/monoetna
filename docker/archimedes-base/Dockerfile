FROM etna-py AS etna-py
FROM development-certs

FROM python:3.8-slim-bullseye
WORKDIR /app

ENV DOCKER_VERSION 19.03.12
ENV PYTHON_VERSION 3.8
ENV DOCKERIZE_VERSION 0.5.0
ENV POETRY_VERSION 1.1.15
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz
ENV PATH "/app/bin:/opt/poetry/bin:$PATH"
ENV POETRY_HOME "/opt/poetry"

COPY --from=1 /certs/rootCA.pem /usr/local/share/ca-certificates/devRoot.crt
RUN chmod 644 /usr/local/share/ca-certificates/devRoot.crt && update-ca-certificates

RUN apt-get update \
   && apt-get install -y \
      curl \
      libopenblas-base \
      libcurl4-openssl-dev \
      libssl-dev \
      wget \
      gnupg \
      bash-completion \
      git \
      tk \
      --no-install-recommends

# add postgres 10 repo for debian bionic
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ bullseye-pgdg main" >> /etc/apt/sources.list.d/pgdg.list
RUN curl https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add -

RUN apt-get update \
   && apt-get install -y \
      postgresql-client-10 \
      --no-install-recommends &&\
    rm -rf /var/lib/apt/lists/*

RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )
RUN curl -sSL https://install.python-poetry.org > /opt/install-poetry.py
RUN python /opt/install-poetry.py --version $POETRY_VERSION

# bash improvements for easier working inside the container
RUN git clone --depth=1 https://github.com/Bash-it/bash-it.git ~/.bash_it && \
    bash ~/.bash_it/install.sh --silent && \
    echo "export SCM_CHECK=false" >> ~/.bashrc

# Add the docker client
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/

# Copy entrypoint scripts
COPY entrypoints /entrypoints
COPY --from=etna-py /etna-py /etna-py
