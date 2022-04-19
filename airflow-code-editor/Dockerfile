FROM node:14
RUN mkdir /airflow-code-editor
COPY . /airflow-code-editor
WORKDIR /airflow-code-editor
RUN npm install && npm run build
