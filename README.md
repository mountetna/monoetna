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

In addition, you'll want to configure some users and projects into your janus environment.

```bash
# Make sure services are already running with make up in top level directory.
make -C janus bash
○ → ./bin/janus add_project 'test-project' 'Test Project'
○ → ./bin/janus add_user some.guy@ucsf.edu Zach Collins password
○ → ./bin/janus permit some.guy@ucsf.edu test-project administrator
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

