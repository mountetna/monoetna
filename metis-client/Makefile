include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

build:
	cargo build --release --target=x86_64-unknown-linux-musl
	cargo install --path .
