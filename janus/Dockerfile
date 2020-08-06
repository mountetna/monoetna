FROM etna-base
ARG APP_NAME
ARG RUN_NPM_INSTALL
ARG SKIP_RUBY_SETUP=1
RUN /entrypoints/build.sh
RUN $(npm bin)/sass --watch lib/client/scss/janus.scss:lib/client/css/janus.css

