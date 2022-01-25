FROM alpine:3.15
RUN apk add --no-cache \
  openssh-client \
  ca-certificates \
  bash
RUN ssh-keyscan github.com > /known_hosts
