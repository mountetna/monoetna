# Monoetna Docker

This directory contains two types of resources:

1. `build_support` -> a set of tools that help build docker images for the environment
2. A number of Dockerfile build directories, generally either for the purpose of supporting common base images
   or providing development related services.

## Images

Monoetna images come in 3 flavors:

1. 'Base images' that are used by other monoetna images.  These base images are either/both shared across multiple
    other images, or used as the image in docker-compose.yml files.
2. 'Application images', typically existing inside a top level directory.  Application images actually get deployed
    and run actively in production.
3. 'Swarm images', typically existing inside the `swarm/` directory.  These images are deployed actively in production,
    but generally are not as actively worked on and may not have tests run in CI to save time.  These images likely
    require manual `PUSH_IMAGES=1 make -C release` invocations to build and release.

## 'Projects'
'projects' are simply directories that exist either at the top level, in docker, or swarm, or in the 
etna/packages directories.  A project can build a docker image *by the same name* simply by having a 
Dockerfile in that directory.  `utils.mk` will know how to build that project automatically when it is discovered
as a dependency in the system.

Dependencies are currently implicitly inferred from two places: docker-compose.yml files and Dockerfiles.
See `utils.mk`#buildable_dependent_image_dockerfiles and buildable_dependent_compose_image_dockerfiles functions.
Basically, in short: __any `FROM` line or `image:` inside of a Dockerfile or docker-compose.yml respectively, constitutes
a project dependency, and `utils.mk` will attempt to find a Dockerfile to build that project's image for you.

Not every project must have a Makefile, but if you want an easy way to execute basic commands, create  `Makefile` and
include atleast `stubs.mk` and likely `docker-compose.mk` to the project to have several useful, project scoped targets
created for you.

## Development vs Production

Development workflows are a bit different than production image workflows, which is important to keep in mind when
understanding how all the docker technologies converge.  A short explanation:

1. Docker images are immutable directory + metadata configurations for applications.  They are built with Dockerfiles
   and are shared both in development and production.
2. docker-compose is a development focused tooling that starts Docker images locally with fake networks and fake
   mounts that simulate a production like environment.  docker-compose is different in production, however, in that
   the source code is not 'bake' and immutable like normal images are.  docker-compose.yml files use *-base images
   and then mount mutable directories containing the mutable code being worked on in development.  This means that
   development is different than production: any local changes you are making, any intermediate files or local installs,
   those will effect your development behavior but not your production behavior.  )
3. docker-swarm: TODO.

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
    top level `Makefile`, several useful development targets that automatically manage the project's docker-compose
   context are created.
    
**IMPORTANT**

By default docker-compose namespaces all its processes, volumes, and networks by the directory name in which
it is invoked.  So that means in general docker-compose will be namespaced with `monoetna_` proceeding it.
If you invoke docker-compose by hand (rather than using the Makefiles, which are recommended), you will need
to override this behavior to ensure your namespace is consistent by either providing `-p monoetna` or `export COMPOSE_PROJECT_NAME=monoetna`.

Recommended that you mainly use the `Makefile` to interface with docker-compose to prevent this from being an issue.
