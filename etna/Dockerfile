FROM etna-base-dev
ENV ETNA_GEM_DEVELOPMENT=1
COPY Gemfile /etna/Gemfile
COPY Gemfile.lock /etna/Gemfile.lock
COPY etna.gemspec /etna/etna.gemspec
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install -j "$(nproc)"

COPY ./package.json /etna/package.json
COPY ./package-lock.json /etna/package-lock.json

# The images tend to build as root, which for host systems is unsafe,
# but in containers is fine.
RUN cd /etna && npm install --unsafe-perm

COPY ./packages/etna-js /etna/packages/etna-js
RUN cd /etna && npm run build


COPY . /etna
