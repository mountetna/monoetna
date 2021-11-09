# monoetna

mono-repository version of etna projects

## Setup for Mac

You'll need `homebrew` and to install a few extra items. Most tools are developed for linux, and thus
there are some minor discrepancies that need to be addressed for mac.

First off, you'll nee to install `homebrew install coreutils` to get the base set of gnu tools.
Secondly, you may need `homebrew install findutils` as well. Try running `type -p gfind || echo need findutils!` to determine
if you need to install findutils as well.

Mostly, the goal is to have standard linux executables.  Macosx unfortunately comes with some out of date BSD variants
of common gnu utils that will cause the command line tools to behave weirdly.

## Directory Structure

Top Level directories fall into one of three types:

1. Deployable applications.
   - `metis` `timur` `edge-apache` and others. Generally these applications fall into one of three categories
     1. `etna` applications. These involve, at minimum, a ruby server, and possibly a javascript asset pipeline. They also share libraries via the `etna` directory.
     2. infrastructure applications. For instance `edge-apache` and `development-psql`.
     3. "Support" applications like `mountetna.github.io` that do not generally need to be run in standard development, and have special properties.
   - Most of these applications define their own `Makefile`, `Dockerfile`, and `docker-compose.yml` files. These
     can provide overwrites of standard behaviors. see `docker/README.md` for more details.
2. Development utilities
   - `bin` contains several useful development scripts for development operations. Some of these are also used by CI.
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

This process can unfortunately be very slow, especially on first build, due to the need to install many many dependencies
between each project and get base images built.  `npm` and `bundle` tend to require heavy, heavy disk and cpu usage
while processing, so expect to let an initial `make up` run for an hour or more.

### Update

Once you've got docker installed and such, you'll want to start by running `make update`.  This command will take a long time,
as it will first need download many things, build our docker images, start databases, update them, install gems and npm,
and a bunch of other filesystem intensive tasks.

Subsequently, after pulling a new set of changes or adjusting Gemfile / package.json / Dockerfiles, run `make update` again
to apply those changes.

`make update` will migrate development and test databases with the current code set.

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
```

### Seeding janus

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
TODO!  A nice way of setting up development timur and magma would be, uh, nice.

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
6. A github action will start and push your changes and deploy to staging.
