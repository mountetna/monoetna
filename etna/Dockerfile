FROM etna-base-dev
ENV ETNA_GEM_DEVELOPMENT=1
COPY Gemfile Gemfile.lock /etna/
RUN cd /etna/ && bundle install
COPY . /etna
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install
