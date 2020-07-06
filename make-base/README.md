# make-base

## What is this?

Provides a 'base' makefile that provides default target implementations and stubs for
projects in order to share development script logic more freely.

Most per-project make files should include one of the following makes based on their needs:

```
# Every make should include.  Ensures a stub target exists for all required projects.
include ../make-base/stubs.mk

# Any project controllable via docker-compose should include
include ../make-base/docker-compose.mk

# Projects that are ruby and use etna should include to receive default sensible testing, migration
# This requires the docker-compose.mk inclusion as well.
include ../make-base/etna-ruby.mk

# Projects that use nodejs should include to receive default sensible testing targets
include ../make-base/node.mk
```

Best documentation is code, however, so check out other projects' `Makefile` for a sense of usage.