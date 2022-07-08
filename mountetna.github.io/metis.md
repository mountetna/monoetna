---
layout: default
group: app
---
# Metis
{:.no_toc}
## File Service
{:.no_toc}

Metis is a file service for Etna applications. It provides the ability to store
binary files in folder hierarchies and access them via HTTP API. The underlying
object storage uses an ordinary filesystem (i.e., files in Metis are stored on
disk as files)

* toc
{:toc}

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

### Browser

Metis provides a browser client that allows file viewing, download, and
upload.

### Command-line

`bin/metis_client` is a command-line client to interact with metis. To run the client requires Ruby 2.5+.

You should install Ruby via [rbenv](https://github.com/rbenv/rbenv) and select one of the above versions (2.5 - 2.6). If you are using a non-Linux operating system, you will need to install the dependencies listed in [the non-Linux OS section](#non-linux-operating-system-dependencies) **before** installing `rbenv` or any of the target Ruby versions.

#### Non-Linux operating system dependencies

When using metis_client on non-Linux operating systems, you may have to manually install some additional libraries before installing `rbenv`. 

  * md5sum
  * readline

For macOS, you can install these both with [homebrew](https://brew.sh/):

```
$ brew install md5sha1sum
$ brew install readline
```

#### Setup

To run the metis client, copy your Janus token (by visiting e.g.
`janus.example.org`) and set this as an environment variable in your shell (e.g.
`export TOKEN=your.janus.token`). Then run `bin/metis_client`.

The first time the client runs it will guide you through setup of these variables:
- metis_host - the host name of the metis server
- metis_uid_name - the token name for your metis_uid
If these values are set incorrectly or not set, the client will fail to connect to Metis.

#### Invocation

With no arguments `metis_client` will put you into an interactive shell where
you may use the commands below.

A valid metis path as optional first argument to `metis_client` will set the
client's project, bucket and folder. E.g., `metis_client
metis://athena/armor/blueprints` will start the shell connected to the 'athena'
project, in the 'armor' bucket, in the 'blueprints' folder.

Subsequent arguments will be treated as commands and executed
non-interactively, e.g. `metis_client metis://athena/armor/blueprints ls` will
list the contents of the `blueprints` folder.

#### Command-line invocation mode

Once you have a set of Metis commands that you want to execute, you can put them into a text file and feed them into metis_client via the command line. This allows you to automate a process with a scheduler, like `cron` or `systemd`.

For example, you might have a text file like the below `instructions.txt`:

```
project my_project
cd data_bucket
get . /local-path-to-copy-data-to
```

Then you can execute these commands against Metis (assuming you have injected the `TOKEN` environment variable)

```
$ metis_client.rb < instructions.txt
```

#### Commands

- *help* - Help! Use this command to get usage information on other commands.
- *project* - List projects or select a project on metis
- *ls* - List a folder path on metis
- *cd* - Change your current metis folder path
- *pwd* - Print your local directory
- *lcd* - Change your local directory
- *get* - Recursively download a directory from metis
- *validate* - Verify that files within a local directory are archived on metis
- *nuke* - Remove files within a local directory that are already archived on metis

If you have edit permission you may also:
- *mv* - Rename a file on metis
- *put* - Upload a folder to metis

#### Avoiding token expiration

Metis downloads may take days, especially for the initial download of a large
corpus. A common issue is token expiration during the download. While
resumption is relatively painless, token expiration will usually interrupt a
workflow. There are some strategies you might employ to avoid token
expiration:

1. metis_client will attempt to juggle your token (keep it current) during an
   active download using Janus's `/refresh_token` endpoint. However, this will
   not re-export a valid token to your shell, so subsequent invocations of
   metis_client (e.g. if you are using command invocation to script downloads)
   will fail.

2. You may use Janus's token generation scheme (see above) to refresh your
   token between downloads. This requires setting a public key on your Janus
   account.

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

## Setup

### Installation

When you update or install the project, make sure to run migrations. You should do this with the built-in command provided in `$ bin/metis migrate` instead of directly using `Sequel`.

To run migrations against your test database, run `$ METIS_ENV=test bin/metis migrate`.

### Data Directory

You can configure where the data blocks are stored in your `config.yml` file. The environment (`test`, `development`, etc.) should have a `:data_path:` value set. Make sure to also create two sub-directories, `uploads` and `data_blocks`. For example, if your `config.yml` looks like:

```yml
:development:
  :data_path: ./data
```

On your disk you'll want to create the following directories:

```sh
data/
  data_blocks/
  named/
  uploads/
```

### Archiving files

When you upload files using the UI, they will appear in your `data_blocks` directory as `temp-` files. To calculate their MD5 hashes (and optionally put them into long-term cloud storage on AWS), you need to run `$ bin/metis archive`.
