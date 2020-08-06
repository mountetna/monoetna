FROM etna-base
WORKDIR /app
ARG APP_NAME
ARG RUN_NPM_INSTALL
ARG SKIP_RUBY_SETUP=1
RUN /entrypoints/build.sh
RUN npm run build
