# monoetna

mono-repository version of etna projects

## Setup for Mac

You'll need `homebrew` and to install a few extra items. Most tools are developed for linux, and thus
there are some minor discrepancies that need to be addressed for mac.

First off, you'll nee to install `homebrew install coreutils` to get the base set of gnu tools.
Secondly, you may need `homebrew install findutils` as well. Try running `type -p gfind || echo need findutils!` to determine
if you need to install findutils as well.

## Directory Structure

Top Level directories fall into one of three types:

1. Deployable applications.
   - `metis` `timur` `edge-apache` and others. Generally these applications fall into one of three categories
     1. `etna` applications. These involve, at minimum, a ruby server, and possibly a javascript asset pipeline. They also share libraries via the `etna` directory.
     2. infrastructure applications. Currently this is mainly the `edge-apache` but could include other vendor applications we intend to build and use in our environment.
     3. "Support" applications like `mountetna.github.io` that do not generally need to be run in standard development, and have special properties.
   - Most of these applications define their own `Makefile`, `Dockerfile`, and `docker-compose.yml` files. These
     can provide overwrites of standard behaviors. see `docker/README.md` for more details.
2. Development utilities
   - `bin` contains several useful development scripts for general operations. Some of these are also used by CI.
   - `development-certs` contains development SSL certs used by the development servers.
3. Build utilities
   - `docker` contains base images and build logic for our docker based build system.
   - `make-base` contains `*.mk` files to be included by deployable applications and provides default `Makefile` command hooks necessary for the build process.
   - `.github` contains yml files defining our github actions based CI processes.

## Docker

Use the `Makefile` at the top of the repo to easily start and stop a set of all etna services that are dockerized.

`make help`
`make up`
`make down`
`make logs`

You will want to use per-project `Makefile`s to access proejct specific databases and bash consoles. eg: `make -C janus help`

### Builds

See more details on how builds behave and how to change the development environment in the `docker/README.md` file.

## Environment setup

You'll want to setup your system's `/etc/hosts` file to support mapping subdomain https urls. The docker
environment runs an apache server that terminates SSL with self signed certs currently. You'll want entries for
each service like the following:

```
127.0.0.1 janus.development.local
127.0.0.1 metis.development.local
127.0.0.1 magma.development.local
127.0.0.1 timur.development.local
127.0.0.1 polyphemus.development.local
127.0.0.1 vulcan.development.local
127.0.0.1 prometheus.development.local
127.0.0.1 airflow.development.local
```

### Seeding janus and metis

In addition, you'll want to configure some users and projects into your janus environment.

A script `./bin/seed_janus` is provided that can setup a default `developer@ucsf.edu` user with pw `password` with
the appropriate privileges to get started, as well as some additional projects.

You can also interact with janus directly to add more projects / users. An example is given below

```bash
# Ensures images are build, bundle has run, and services come up
make up
# Enters bash console in janus
make -C janus bash
○ → ./bin/janus add_project 'test-project' 'Test Project'
○ → ./bin/janus add_user developer@ucsf.edu Developer LastName password
○ → ./bin/janus permit developer@ucsf.edu test-project administrator
○ → ./bin/janus add_project ipi 'Immuno Profiler Project'
○ → ./bin/janus permit developer@ucsf.edu ipi administrator
```

### Seeding timur and magma

You'll want to also likely setup some useful data for timur and magma. These seeds are large (~2GB total) but give you
fairly useful data to test again.

```bash
# Make sure other services are not accessing the db by turning them off
make down
export TOKEN=JANUSTOKEN
./bin/seed_databases
```

NOTE: This does delete your development timur and magma databases on each run. It is a complete in place
replacement of those with the seed data.

## Debugging issues

First off, check that your processes started correctly and are running.

```
make ps
```

Logs can be useful. An aggregate of all project logs can be seen with

```
make logs
```

But this may be overwhelming. Try narrowing it down to a specific project with

```
make -C timur logs
```

Is the problem something to do with docker, rather than the running processes?

Because docker's namespace is global to the host, it's possible for networks and volumes it creates to conflict.
If you need to, try removing a network or volume to resolve these conflicts. They will be auto re-created on the next `make up`

```
docker volumes
docker network ls
docker volume rm NAME
docker network rm NAME
```

## Workflow

1. Create a branch in monoetna as per your feature.
2. Write code that cross cuts all repositories
3. Push your changes to master, aka the monoetna repo.
4. Create PR in monoetna repo, showing the entire change set
5. Merge to master on monoetna
6. A github action will start and push your changes to the subtree repos automatically!

### Branches

Both `./bin/pull-subtrees` and `./bin/push-subtrees` can push to branches
that exist on child repositories. However, in that case, _only repositories
that actually already contain a branch matching the name of the current
monoetna branch will be pushed or pulled_.

IE:
etna has a branch named zc/my-feature
metis has a branch named zc/my-feature
etna-js does not

If I `./bin/push-subtrees` from branch zc/my-feature, only etna and metis will have changes pushed to.

In general, it is wiser to focus on master level development.
