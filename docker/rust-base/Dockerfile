FROM rust:latest
RUN apt-get update
RUN apt-get install musl-tools -y
RUN rustup target add x86_64-unknown-linux-musl
ENV RUSTFLAGS="-Clinker=musl-gcc"

WORKDIR /app

COPY Cargo.lock .
COPY Cargo.toml .
COPY src src

RUN cargo build --target=x86_64-unknown-linux-musl
RUN cargo build --release --target=x86_64-unknown-linux-musl
RUN rm -rf ./target/x86_64-unknown-linux-musl/release/deps/base-deps*
RUN rm -rf ./target/x86_64-unknown-linux-musl/debug/deps/base-deps*
