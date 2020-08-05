# monoetna
mono-repository version of etna projects

## Docker

Use the `Makefile` at the top of the repo to easily start and stop a set of all etna services that are dockerized.

`make help`
`make up`
`make down`
`make logs`

You will want to use per-project `Makefile`s to access proejct specific databases and bash consoles. eg: `make -C janus help`

### Docker environment setup

You'll want to setup your system's `/etc/hosts` file to support mapping subdomain https urls.  The docker
environment runs an apache server that terminates SSL with self signed certs currently.  You'll want entries for
each service like the following:

```
janus.development.local 127.0.0.1
metis.development.local 127.0.0.1
magma.development.local 127.0.0.1
timur.development.local 127.0.0.1
```

### Seeding janus and metis

In addition, you'll want to configure some users and projects into your janus environment.

```bash
# Ensures bundle has run
make prepare
# Enters bash console in janus
make -C janus bash
○ → ./bin/janus add_project 'test-project' 'Test Project'
○ → ./bin/janus add_user developer@ucsf.edu 'Developer' password
○ → ./bin/janus permit developer@ucsf.edu test-project administrator
```

To use magma and timur effectively, you'll want to create an `ipi` project as well.

```bash
# Make sure services are already running with make up in top level directory.
make -C janus bash
○ → ./bin/janus add_project ipi 'Immuno Profiler Project'
○ → ./bin/janus permit developer@ucsf.edu ipi administrator
```

### Seeding timur and magma

You'll want to also likely setup some useful data for timur and magma.  These seeds are large (~2GB total) but give you
fairly useful data to test again.

```bash
# Make sure other services are not accessing the db by turning them off
make down
./bin/seed_databases
```

NOTE:  This does delete your development timur and magma databases on each run.  It is a complete in place
replacement of those with the seed data.

## Issues

If you run into issues building docker images, one quick solution is to try

```bash
cd docker
docker-compose -p monoetna --no-cache
```

Also, check if your processes are having issues with

```bash
make ps
make logs
make -C metis logs
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
that exist on child repositories.  However, in that case, *only repositories
that actually already contain a branch matching the name of the current
monoetna branch will be pushed or pulled*.

IE:
etna has a branch named zc/my-feature
metis has a branch named zc/my-feature
etna-js does not

If I `./bin/push-subtrees` from branch zc/my-feature, only etna and metis will have changes pushed to.

In general, it is wiser to focus on master level development.

### Future

To simplify the workflow, we could configure github actions to be capable
of pushing and pulling automatically.

