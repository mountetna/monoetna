# Monoetna Docker

This directory contains the nuts and bolts necessary to build monoetna's docker images for local
development and remote deployment.


## Images

Monoetna images come in 3 flavors:

1. 'Base images' that are used by other monoetna images.  These belong in the `docker/` directory as subdirectories.
2. 'Application images', which are simply top level directories that contain a `Dockerfile` and appropriate `Makefil`.  These often use a base image.
2.  Remote images that are downloaded from dockerhub.   They simply fall out of the above two cases.

### Creating a new Base Image

Create a new sub directory in `docker`, and add a `Makefile` to build that image.
Copying from one of the other available images is a good starting place.  Use `build_image` script for consistency.

### Creating a new application image

Create a directory with a unique name in the top level of the monoetna project.
Add a new `Makefile` with contents like so:

```
include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
include ../make-base/etna-ruby.mk
include ../make-base/node.mk
```

`make-base/README.md` has more details on how to build our your makefile.

From there, you can use `make -C config-ready` to generate a base Dockerfile and docker-compose.yml file with
further options to configure.

## docker-compose.yml

There are a few variants of docker-compose.yml files that are used exlusively to manage the development
processes.

1. `monoetna/docker-compose.yml` Is simply a link to `monoetna/docker/docker-compose.yml` allowing for convenient 
    invocation of docker-compose commands from the root directory.
2. `monoetna/docker/docker-compose.yml` is generated in `docker/Makefile` anytime an application's
    `docker-compose.yml` is modified.  ESsentially, it uses the `docker/compose` script to join together each
    individual projects' compose.yml files together into one large multi-process one.
3. Application level `docker-compose.yml` files.  By including `make-base/docker-compose.mk` in a project's
    top level `Makefile`, this file becomes managed by the `docker/default_compose` script, which provides some
    defaults and kept up to date via merging logic with common mixin files.
4. Adhoc `docker-compose.yml` files.  These may represent applications that do not make use of `make-base/docker-compose.mk`
    and instead configure manually their settings, but are still included in the top level compose file.
    
**IMPORTANT**

By default docker-compose namespaces all its processes, volumes, and networks by the directory name in which
it is invoked.  So that means in general docker-compose will be namespaced with `monoetna_` proceeding it.
If you invoke docker-compose by hand (rather than using the Makefiles, which are recommended), you will need
to override this behavior to ensure your namespace is consistent by either providing `-p monoetna` or `export COMPOSE_PROJECT_NAME=monoetna`.

Recommended that you mainly use the `Makefile` to interface with docker-compose to prevent this from being an issue.
