# monoetna

mono-repository version of etna projects

## Setup for Mac

### Getting started 

1. Build all the images! Run: `make -f Makefile-mac build-dev-etna-images` and `make -f Makefile-Mac build-airflow-images`
2. Spin up the etna container and install ruby/js dependencies  `make -f Makefile-mac etna-libs-ruby` and `make -f Makefile-mac etna-libs-js`.
This installs JS dependencies locally at `etna/node_modules/` and then this directory is subsequently mounted into containers.
It is unclear where the gems are installed... (TODO: look into this)
3. Spin up the webapps: Run: `make -f Makefile-mac web-up` and `make -f Makefile-mac airflow-up`
4. Run migrations for the webapps `make -f Makefile-mac migrate-all`


### Hardcoded values you must change

On line 23 of `etna/docker-compose.yml we have:

`/home/home/etna.yml:/root/etna.yml`

You must change the first part of the path mapping to:

`/my-home-directory/etna.yml:/root/etna.yml`

The etna gem assumes that the etna.yml file exists in your home directory.

### Example project

First you must install the etna gem locally, and then create a etna.yml file.

#### Creation

1. Create some example projects in janus `make -f Makefile-mac janus-seed` # Is this needed?
2. Create the example projects in magma `make -f Makefile magma-create-project`

#### Seeding some models

If you would like to download the `example` project models in production and use them locally:

1. Copy your PROD janus token `make -f Makefile-mac magma-copy-example-models`
2. Copy your dev janus token `make -f Makefile-mac magma-create-example-models`

#### Seeding some data

TODO

#### Vulcan

The above instructions do not build Vulcan. To build vulcan you must build archimedes images: `make -f Makefile-mac build-dev-archimedes-images`.
This takes a really long time

## Directory Structure

Top Level directories fall into one of four types:

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
4. Swarm Applications
   - Applications under the `swarm` directory are deployed into infrastructure.  In most ways, they behave like top level directory apps except in one crucial detail:
these apps do not run their tests in CI during a branch PR.   The intention is that these apps do not change as often as other apps, thus spending CI resources
building them for every PR is not ideal.  That said, `release-test` is still run for these applications in `staging` and `production` during generation of assets.
The idea here is that branch build times should be slightly faster since they iterate and change more often than mainline staging and production.
If you make changes to these applications, you'll want to test and verify locally.  See *Build System* below.

## Build System

### Linux

The monoetna repo comprises multiple modular applications and libraries that may have one or even more than one different
final locations.  The build system helps support a consistent interface for generating and deploying our system in this
modular way using a singular repository.

Generally, the easiest way to interact with the build tools is through `make` in the repository's top level directory,
although each project contains its own make file and can be interacted with in its own directory as well.

The top level `Makefile` contains several general high level targets, most notably:

```
make up # starts the docker-compose.yml processes containing ALL development apps
make down # stops any and all docker-compose.yml processes associated with this repo
make update # runs local migrations, gem, npm, poetry installs, and any other local development maintenance tasks.
make ps # lists the docker-compose processes running.
```

Some tasks are better oriented for running in the context of specific apps.  you can use `make -C magma` for example to
run commands specifically in a directory (application) context.  A few notable commands useful in application contexts are

```
make -C metis bash # opens a bash process in the docker-compose context of metis app
make -C metis start # starts ONLY the metis processes in docker-compose
make -C metis release-test # runs all tests for metis
make -C metis release-build # builds all images for metis
PUSH_IMAGES=1 IMAGE_PREFIX=etnaagent/ IMAGE_POSTFIX=:production make -C metis release # builds and releases metis to production!  Use only for hot deployments.
```

Greater detail for the implementation and available commands exists in the `make-base` directory.

### Builds

When a branch is pushed to github, CI will execute the contents of `.github/branches.yml`.  Generally, this will execute
all tests for most (but not all) projects, validating that your changeset is safe to merge.

When your branch is merged into master, `.github/master.yml` will execute a `make release` in the top directory, *thereby
building EVERY project and pushing the resulting artifacts to dockerhub*.  However, master.yml tags these images with `:staging`,
so they only effect services which are running in staging mode.

After merging into the `production` branch (usually via a PR from `staging`, but in an emergency you can force push anything to `production` branch),
`.github/production.yml` will kick off a `make release` that pushes to dockerhub with the tag `:production`.

Services in production are constantly monitoring for image updates in dockerhub matching their tags.  When they find one,
they will generally automatically pull and update themselves from that new image.

Generally.  But not all services will do this.  Check the service definitions in`swarm/`.  Services whose labels include
`autoupdate=true` are those that will pull and update automatically from github.  Otherwise, you will need to use the
portainer ui and force an image pull update there.

## Environment setup

When runnin the docker-compose development environment, an apache instance is launched with self signed certs at specific
local host names given by the configuration in `development-edge-httpd.conf`.  However, for your browser to actually accept
and resolve these local hostnames, you should configure your system's dns host resolution to accept them.

This can vary by system, but in general, linux like systems should allow direct entry to the `/etc/hosts` file.  You'll
want to maintain entries like the following based on which services you are testing:

```
127.0.0.1 janus.development.local
127.0.0.1 metis.development.local
127.0.0.1 magma.development.local
127.0.0.1 timur.development.local
127.0.0.1 polyphemus.development.local
127.0.0.1 vulcan.development.local
127.0.0.1 prometheus.development.local
```

### Gemfiles

Currently, all etna ruby projects derive from the `etna` base image, which shares a single `Gemfile` and `Gemfile.lock`.
This is intentional design for performance, storage, and maintenance simplicity.

Adding new gems, one should modify the `etna/Gemfile`, then update the lock file with:

```bash
make -C etna bash
> bundle install
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

### Duplicating magma schemas

You can copy and apply magma project modeling to your local environment to support testing against production like schemas.
Firstly, you'll need a properly configured `etna.yml` file for both the production environment and your local environment.
Next, you'll need the `etna` gem installed on a local ruby.  You can use `bin/reinstall-gem` to do this from your local sources.

```
# export your production token
export TOKEN=asjlkjadsf

# Copies production mvir1 into a csv file
$(gem env path)/bin/etna administrate models copy_template mvir1 --file mvir1_template.csv --environment production

# Swap out for your development janus token here!
export TOKEN=mkaljerljdas
$(gem env path)/bin/etna administrate models apply_template mvir1 --file mvir1_template.csv --environment development
```

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
