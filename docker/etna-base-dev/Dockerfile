FROM development-certs

FROM ruby:2.7.6-bullseye
WORKDIR /app

ENV DOCKER_VERSION 20.10.6
ENV RUBY_VERSION 2.7.6
ENV RUBY_MAJOR_VERSION 2.7.0
# This is so dumb -- but if BUNDLER_VERSION is set as an environment variable,
# it causes issues with bundler... bleck.  So, so add an underscore at the end.
ENV BUNDLER_VERSION_ 2.3.15
ENV DOCKERIZE_VERSION 0.5.0
ENV NODE_VERSION 14.17.1
ENV DOCKERIZE_URL https://github.com/jwilder/dockerize/releases/download/v${DOCKERIZE_VERSION}/dockerize-linux-amd64-v${DOCKERIZE_VERSION}.tar.gz
ENV RACK_ENV development

ENV BUNDLE_PATH="/bundle"
ENV BUNDLE_BIN="$BUNDLE_PATH/ruby/$RUBY_MAJOR_VERSION/bin"
ENV BUNDLE_SILENCE_ROOT_WARNING=0
ENV BUNDLE_APP_CONFIG="/etna/.bundle"
ENV BUNDLE_GEMFILE="/etna/Gemfile"
ENV PATH "/app/bin:/etna/node_modules/.bin:/app/node_modules/.bin:$BUNDLE_BIN:$PATH"

COPY --from=0 /certs/rootCA.pem /usr/local/share/ca-certificates/devRoot.crt
RUN chmod 644 /usr/local/share/ca-certificates/devRoot.crt && update-ca-certificates

# add postgres 10 repo for debian bionic
RUN echo "deb http://apt.postgresql.org/pub/repos/apt/ bullseye-pgdg main" >> /etc/apt/sources.list.d/pgdg.list
RUN curl https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add -

RUN apt-get update \
   && apt-get install -y \
      curl \
      libopenblas-base \
      libcurl4-openssl-dev \
      libssl-dev \
      parallel \
      postgresql-client-10 \
      rsync \
      sshpass \
      libvips \
      libvips-tools \
      sqlite3 \
      --no-install-recommends &&\
    rm -rf /var/lib/apt/lists/*

RUN curl https://nodejs.org/dist/v${NODE_VERSION}/node-v${NODE_VERSION}-linux-x64.tar.gz | tar --strip-components=1 -xzC /usr/local/ \
    && npm install -g yarn typescript

RUN curl -o /tmp/dockerize.tgz -L $DOCKERIZE_URL && ( cd /usr/bin && tar xzf /tmp/dockerize.tgz )

RUN rm /usr/local/lib/ruby/gems/2.7.0/specifications/default/bundler-*.gemspec
RUN gem install bundler --default -v "=$BUNDLER_VERSION_"
## !! Note the rubygems version locked down here is meant to just be a temporary fix and should be removed when we do update Ruby !!
RUN gem update --system 3.4.22

# We also the GEM_PATH to be this BUNDLE_PATH so that ides (Rubymine) will find gems easier.
ENV GEM_HOME="$BUNDLE_PATH"
ENV GEM_PATH="/root/.gem/ruby/$RUBY_MAJOR_VERSION:/usr/local/lib/ruby/gems/$RUBY_MAJOR_VERSION:/usr/local/bundle:$BUNDLE_PATH/ruby/$RUBY_MAJOR_VERSION"

# bash improvements for easier working inside the container
RUN git clone --depth=1 https://github.com/Bash-it/bash-it.git ~/.bash_it && \
    bash ~/.bash_it/install.sh --silent && \
    echo "export SCM_CHECK=false" >> ~/.bashrc

# Add the docker client
RUN curl https://download.docker.com/linux/static/stable/x86_64/docker-$DOCKER_VERSION.tgz > /tmp/docker.tgz
RUN tar -xf /tmp/docker.tgz -C /tmp/
RUN mv /tmp/docker/* /usr/local/bin/

RUN mkdir -p /opt
COPY entrypoints /entrypoints

#RUN bundle install

COPY bin/* /bin/
RUN chmod +x /bin/*
