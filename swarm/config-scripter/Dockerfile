FROM node:17-slim
RUN mkdir /app
WORKDIR /app
COPY package-lock.json /app/
COPY package.json /app/
RUN npm install
RUN ls -la /