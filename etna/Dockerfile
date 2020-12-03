FROM etna-base
ENV ETNA_GEM_DEVELOPMENT=1
RUN rm -rf /app && ln -sfT /etna /app
RUN bundle install

# Install R from CRAN repo, to get latest
RUN gpg --keyserver keys.gnupg.net --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add -
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update && \
    apt-get install -y r-base \
    r-base-dev
