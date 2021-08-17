FROM python:3.8-slim
WORKDIR /app

ENV DOCKER_VERSION 19.03.12
ENV PYTHON_VERSION 3.8
ENV DOCKERIZE_VERSION 0.5.0
ENV POETRY_VERSION 1.1.4
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz
ENV PATH "/app/bin:/opt/poetry/bin:$PATH"
ENV POETRY_HOME "/opt/poetry"

COPY certs/rootCA.pem /usr/local/share/ca-certificates/devRoot.crt
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
      libicu63 \
      libreadline7 \
      --no-install-recommends

# add postgres 10 repo for debian stretch
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ stretch-pgdg main" >> /etc/apt/sources.list.d/pgdg.list
RUN wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add -

RUN apt-get update \
   && apt-get install -y \
      postgresql-client-10 \
      --no-install-recommends &&\
    rm -rf /var/lib/apt/lists/*

RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py > /opt/install-poetry.py
RUN python /opt/install-poetry.py --version $POETRY_VERSION

# Install R from CRAN repo, to get latest
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '0xE19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update && \
    apt-get install -y r-base \
    r-base-dev

# Install dittoSeq into R
# This is a long process, so bake it in the image
#       instead of using the entrypoints.
# Taking out for performance. Leaving here as comment in case
#       we decide we need dittoSeq back.
# RUN echo 'if (!requireNamespace("remotes", quietly = TRUE))\n    install.packages("remotes")\n\nremotes::install_github("dtm2451/dittoSeq")' > install_dittoseq.txt
# RUN R < install_dittoseq.txt --save

# bash improvements for easier working inside the container
RUN git clone --depth=1 https://github.com/Bash-it/bash-it.git ~/.bash_it && \
    bash ~/.bash_it/install.sh --silent && \
    echo "export SCM_CHECK=false" >> ~/.bashrc

# Add the docker client
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/

# Copy entrypoint scripts
COPY src/entrypoints /entrypoints
COPY magby /magby

