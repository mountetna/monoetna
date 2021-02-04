FROM archimedes-base
# Perform these steps first to allow better caching behavior
#RUN bundle install
WORKDIR /app
ADD src/poetry.lock /app/poetry.lock
ADD src/pyproject.toml /app/pyproject.toml
RUN poetry install
ADD src /app/
RUN /entrypoints/build.sh
