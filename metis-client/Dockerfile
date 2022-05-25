FROM rust-base
COPY Cargo.lock .
COPY Cargo.toml .
COPY src src
RUN cargo build --release --target=x86_64-unknown-linux-musl
RUN cargo install --path .