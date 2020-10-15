FROM etna-base
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install
