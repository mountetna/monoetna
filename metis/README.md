![Run tests](https://github.com/mountetna/metis/workflows/Run%20tests/badge.svg)

# Metis

Metis is a file service for Etna applications. It provides the ability to store
binary files in folder hierarchies and access them via HTTP API. The underlying
object storage uses an ordinary filesystem (i.e., files in Metis are stored on
disk as files)

##

### Docker

#### Usage

A `Makefile` is present to simplify usage of docker via docker-compose.
You will need to ensure that
1. `docker` is installed and running its daemon. `docker ps` should succeed.  (Your user may need to belong to the `docker` group)
2. `docker-compose` is installed and available on your path.  `docker-compose --help` should succeed.
3.  Port 3000 is available on your host machine to bind to IFF you are accessing it via `https://localhost:3000`

`Makefile` is used over other task managers because virtually every modern system can get `gnumake` atleast fairly easily,
without having to install other software outside the containers.

`make help` lists all commands for using the docker system.
* `make up` brings up the database, webpack, and ruby server
* `make down` brings it all down
* `make migrate` runs migrations against the database
* `make logs` shows logs of all running containers
* `make bash` opens a bash shell in the context of the ruby server
* `make psql` connects a to the metis container database

#### Files

The main files related to running metis in development are:

* config.yml(.template) - Contains runtime configuration for your environment.  
The config.yml.template provides a default set of configuration that will work
with the docker setup.  On first `make up` it will be copied to `config.yml`.  Note
that `config.yml` itself is not checked in, so you can configure as you need locally.
* docker-compose.yml - defines the services that will run together
for the metis project specifically.  Edit this to modify any environment
variables, or docker setup, specific to metis's development environment.
* docker/app/Dockerfile - a dockerfile that sets up and runs the metis web app.  Configures
ruby, installs bundler, loads all files from the docker/app/* into the container.
* docker/app/docker-entrypoint.sh - A script that is run before any command inside
the container.  Ensures gems are loaded up to date, waits for the db to be available,
performs potential npm installation, and then runs the command.  Used both by the
webserver and the webpack bundle.
* docker/app/puma.sh - the development puma server command.
* docker/db/Dockerfile - The relatively simple docker file that brings up the postgres db
for metis.
* docker/db/docker-entrypoint-initdb.d - a directory containing scripts that will be
run to initialize the database container, mostly just initializes the actual databases
and creates roles.


## Organization

### Projects

As with all Etna applications, the basic organizational unit of Metis is the
project - the user only has rights to data from a project according to their
project role.

### Buckets

Buckets are the root-level containers for each project. A bucket is access
restricted, either by role or by access list. Bucket names are restricted to
`symbol_names`, i.e. `[A-Za-z0-9_]+`.

### Folders

Folders are collections of files and folders. Folders may be "protected",
preventing the creation of new files or folders.

### Files

Files may have any content. Metis will compute MD5 sums for each file and (if
available) back them up using cloud storage. Files may also be protected from
modification by admins.

### Names

File (and folder) names must match `[^<>:;,?"*\|\/\x00-\x1f]+`, i.e. excluding
common wildcard, separator and control characters.

## Client

Metis provides a browser client that allows file viewing, download, and
upload.

## API

### Listing

Listing a folder path will return a JSON list of files and folders at that
path, including an HMAC-signed `download_url`:

```
GET /:project_name/list/:bucket_name/*folder_path
```

### Buckets

You may get a list of visible buckets for your project:
```
GET /:project_name/list/
```

Admins may create, update or delete a bucket:
```
POST /:project_name/bucket/create/:bucket_name { owner, access, description }
POST /:project_name/bucket/update/:bucket_name { access, description, new_bucket_name }
DELETE /:project_name/bucket/remove/:bucket_name
```

Bucket access is either a role `administrator, editor, viewer` or a
comma-separated list of Etna user ids (emails).

### Folders

Editors may create, remove and rename folders.
```
POST /:project_name/folder/create/:bucket_name/*folder_path
DELETE /:project_name/folder/remove/:bucket_name/*folder_path
POST /:project_name/folder/rename/:bucket_name/*folder_path { new_folder_path }
```

Admins may protect or unprotect folders:
```
POST /:project_name/folder/protect/:bucket_name/*folder_path
POST '/:project_name/folder/unprotect/:bucket_name/*folder_path
```

### Files

Editors may remove or rename a file:
```
DELETE /:project_name/file/remove/:bucket_name/*file_path
POST /:project_name/file/rename/:bucket_name/*file_path { new_file_path }
```

Admins may protect or unprotect a file:
```
POST /:project_name/file/protect/:bucket_name/*file_path
POST /:project_name/file/unprotect/:bucket_name/*file_path
```

### Uploads

The basic upload cycle first requires the editor to authorize a new file upload
into a project, bucket and path with their Etna auth token:
```
POST /authorize/upload { project_name, bucket_name, file_path }
```
If the intended file path is invalid, or the destination is locked, the
upload authorization will fail.

This returns an HMAC-signed URL to the upload endpoint, which may be
used to perform the upload:
```
POST /:project_name/upload/:bucket_name/*file_path { action, ... }
```

This endpoint can perform three actions:

We initiate the upload by communicating what we intend to send:
```
{ action: 'start', file_size, next_blob_size, next_blob_hash }
```

We repeatedly send blobs of binary data (using multipart post):
```
{ action: 'blob', blob_data, next_blob_size, next_blob_hash }
```

If the data does not hash correctly upon receipt, Metis will reject the blob.

When we have sent the final blob, the upload completes and Metis returns the
newly-minted file JSON.

If we are dissatisfied with our progress we may cancel the upload:
```
{ action: 'cancel' }
```

## Documentation

For the most current documentation about Metis, please refer to the [Mount Etna documentation blog](https://mountetna.github.io/).