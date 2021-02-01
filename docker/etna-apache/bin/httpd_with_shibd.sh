#!/usr/bin/env bash

set -e

/etc/init.d/shibd start
httpd -DFOREGROUND

