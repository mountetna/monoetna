# How the pipeline works

## State management

- volume that contains files_to_update.txt
- runs.state 

## Pipeline details

### Job 1 - SFTP File discovery

- checks runs.state.last_scan in the db 
- if it is not set, we set start_time to the initial_start_scan_time
- if it is set, we set start_time to the last_scan timestamp
- if the config.interval param is set, we set end_time to the last_scan timestamp + interval
- if the config.interval param is not set, we set end_time to the current time
- performs a regex search on the ftp server for all files that > start_time and < end_time
- grab their file name and timestamp it which it was updated and write it to a file:
called:  {file_path}/{argo_id}/files_to_update.txt - in csv form.
- writes number of files_to_update to the db

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