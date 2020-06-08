![Run tests](https://github.com/mountetna/metis/workflows/Run%20tests/badge.svg)

# Metis

Metis is a file service for Etna applications. It provides the ability to store
binary files in folder hierarchies and access them via HTTP API. The underlying
object storage uses an ordinary filesystem (i.e., files in Metis are stored on
disk as files)

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