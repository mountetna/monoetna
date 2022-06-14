FROM etna-base-dev
ENV ETNA_GEM_DEVELOPMENT=1
COPY . /etna
RUN rm -rf /app && ln -sfT /etna /app

RUN bundle config set --local no_prune 'true'
RUN bundle config set --local deployment 'true'
RUN bundle install -j "$(nproc)"
RUN cd /etna && npm run build