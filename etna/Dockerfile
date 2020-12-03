FROM etna-base
ENV ETNA_GEM_DEVELOPMENT=1
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install

# Install R from CRAN repo, to get latest
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update && \
    apt-get install -y r-base \
    r-base-dev
