# How the pipeline works

## State management

- volume that contains files_to_update.txt
- db table `cat_ingestion`: 
	- id
	- argo_id
	- last_scan
	- num_files_to_update
	- num_c4_files_updated
	- num_metis_files_updates
	- updated_at
	- modified_at

## Pipeline details

### Job 1 - File discovery

- checks cat_ingestion.last_scan in the db 
- performs a regex search on the ftp server for all files that > last_run
- grab their file name and md5 and write it to a file files_to_update.txt - in csv form file, hash
- writes number of files_to_update 

### Job 2 - C4 update

- wakes up and reads the files_to_update.txt
- connects to the c4 and searches for existing files in the filesystem from the csv
- if they exists update them (download from ftp server and upload)
- for all new files just write them (download from ftp server and upload)
- write num_c4_files_updated to the db
- if any files fail to update, write a file: c4_failed_files.txt, return error

### Job 3 - Metis update

- wakes up and reads the files_to_update.txt
- queries metis api and searches for existing files from the csv
- if they exists update them (download from ftp server and upload)
- for all new files just write them (download from ftp server and upload)
- write num_metis_files_updated to the db
- if any files fail to update, write a file: metis_failed_files.txt, return error

### Job 4 - Reporting

- once job 3 and 2 are finished:
	- compute
		num_files_to_update -  num_c4_files_updated
		num_files_to_update -  num_metis_files_updated
	- report the above to slack

#### Parallelism

Job 2 and Job 3 can run in parralel