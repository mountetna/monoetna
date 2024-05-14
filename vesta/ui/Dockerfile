# TODO: consider switching to next.js-provided Dockerfile
# https://github.com/vercel/next.js/blob/canary/examples/with-docker/Dockerfile
FROM node:20.11

WORKDIR /app

COPY package.json package-lock.json ./
RUN npm ci

COPY . .
RUN npm run build

ENV NEXT_TELEMETRY_DISABLED 1
ENV NODE_ENV production

EXPOSE 3000
ENV PORT 3000
ENV HOSTNAME "0.0.0.0"

CMD ["npm", "start"]