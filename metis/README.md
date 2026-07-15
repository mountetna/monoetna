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
POST /:project_name/folder/unprotect/:bucket_name/*folder_path
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

## File Lifecycle

Every file in Metis is backed by a **DataBlock** — the physical content on disk — while the file record is a named pointer to that block. The two can evolve independently, which is what makes deduplication and safe deletion possible.

### 1. Upload completes

When an upload finishes, a DataBlock is created with a temporary identifier since the real checksum is not yet known. The uploaded file is moved to a holding location on disk and the file record is immediately associated with this temporary block. A `create_datablock` and a `link_file_to_datablock` event are logged.

### 2. Checksum computed

A background checksum job runs after the upload and computes the real MD5 of the file. One of three outcomes follows:

- **New content** — the temporary block is promoted to a permanent one with its real checksum and moved to its final location on disk. A `resolve_datablock` event is logged.

- **Duplicate of existing live content** — the file is re-pointed to the already-existing block. The temporary block is discarded. A `reuse_datablock` event is logged. This is deduplication.

- **Duplicate of previously vacuumed content** — the file data is moved back to the old block's location, the old block is brought back to life, and the file is re-pointed to it. The temporary block is discarded. A `restore_datablock` event is logged. This is resurrection — the original block's history is preserved and it becomes active again.

### 3. File copied or moved

When a file is copied to a new path, or linked to a new location, a `link_file_to_datablock` event is logged. Multiple files can point to the same DataBlock at the same time.

### 4. File deleted

When a file is removed, the file record is deleted and an `unlink_file_from_datablock` event is logged. The underlying DataBlock is not touched — it may still be referenced by files in other projects.

### 5. DataBlock vacuumed

Once a DataBlock is orphaned and eligible for vacuum, its physical file is deleted from disk and the block is marked as removed. Importantly, the block's database record is never destroyed — this is what makes resurrection possible if the same content is uploaded again later. A `remove_datablock` event is logged.

### 6. DataBlock blacklisted

A DataBlock can be marked as blacklisted to signal that its content was intentionally purged, distinct from a normal vacuum. Blacklisted blocks are excluded from active use but their records are retained for audit purposes.

---

## Vacuum Eligibility

A DataBlock is eligible for vacuum only when all of the following are true:

1. **It was linked to a file at some point** — the ledger has a record of the block being associated with a file in this project, whether through a fresh upload, deduplication, or restoration.

2. **It was later unlinked** — the ledger has a record of the block's last file being removed in this project.

3. **No file anywhere currently points to it** — not just in this project, but in any project. If another project still has a live file backed by this block, it cannot be vacuumed.

4. **It has not already been vacuumed since its last restoration** — a block that has been vacuumed is normally excluded from future vacuum runs. However, if it was subsequently restored (because the same content was re-uploaded), that exclusion is lifted and it becomes eligible again.

> **Temp blocks are never vacuumed.** A block that is still in the middle of checksum resolution is never eligible for vacuum, preventing a race condition where an in-flight upload's data could be deleted before the checksum job finishes.

### Blocked datablocks

A DataBlock can be **orphaned by one project but still live in another**. For example, if project A and project B both have files backed by the same block, and project A deletes its file, the block is orphaned from project A's perspective but project B still needs it. In this case the block is considered **blocked** — it cannot be vacuumed until all projects have orphaned it.

The stats endpoint (see below) surfaces which of your orphaned blocks are blocked and which projects are holding them, so you can coordinate with those teams before vacuuming.

---

## Ledger

The ledger records all data block events in Metis and is the source of truth for vacuum decisions. There are two modes of operation — these are not controlled by a flag, they are just useful conceptual distinctions.

### Backfilled

Backfilled mode refers to datablocks and events that existed before the ledger was enabled system-wide. These were backfilled using the `SYSTEM_BACKFILL` trigger. Only link and unlink events can be reliably reconstructed during backfilling.

Two commands must be run to fully backfill a project:

```
bin/metis backfill_data_block_ledger --project_name <project_name> --links
```

Backfills link events for a given project. Run this for every project.

```
bin/metis backfill_data_block_ledger --orphaned
```

Finds all orphaned datablocks across all projects and creates unlink events for them. No project name is needed because datablocks can span multiple projects and the original project is often unknown.

It is best to run `--links` before `--orphaned`. Vacuuming backfilled datablocks only requires the unlink events created by `--orphaned`, but running `--links` is logically nice and is necessary for future tracked-mode vacuums.

### Tracked

Tracked mode is the standard mode of operation once the ledger is enabled system-wide and backfilling is complete. Events are recorded in real-time as files are created, linked, unlinked, and so on.

### Commands

#### `bin/metis ledger_stats`

Shows ledger event counts and vacuum candidates for a project. Run this before vacuuming to understand what is ready and what is blocked.

```
bin/metis ledger_stats --project_name <project_name>   # tracked mode
bin/metis ledger_stats --backfilled                    # backfilled mode
```

Output includes how many datablocks are ready to vacuum, how many are blocked and by which projects, and a detailed list of the orphaned datablocks with their associated file paths.

#### `bin/metis vacuum_datablocks`

Deletes orphaned datablocks. Defaults to a dry run — add `--commit` to actually delete. The CLI will prompt for confirmation before committing.

```
bin/metis vacuum_datablocks <project_name>          # dry run
bin/metis vacuum_datablocks <project_name> --commit # actually delete
bin/metis vacuum_datablocks backfilled --commit     # vacuum backfilled orphans
```

#### `bin/metis backfill_data_block_ledger`

Reconstructs ledger history for datablocks that existed before the ledger was enabled. See the Backfilled section above for the full backfilling workflow.

## Documentation

For the most current documentation about Metis, please refer to the [Mount Etna documentation blog](https://mountetna.github.io/).