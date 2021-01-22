FROM vulcan-base
# Perform these steps first to allow better caching behavior
#RUN bundle install
ADD src /app/
ARG APP_NAME
RUN /entrypoints/build.sh
