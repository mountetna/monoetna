FROM archimedes-base
# Perform these steps first to allow better caching behavior
#RUN bundle install
WORKDIR /app
COPY src/poetry.lock /app/poetry.lock
COPY src/pyproject.toml /app/pyproject.toml
RUN poetry install
ARG FULL_BUILD=1
COPY src /app/
RUN /entrypoints/build.sh
