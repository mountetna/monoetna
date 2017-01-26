## Upload Cycle Description.

I wanted to add a note here to describe the full upload cycle. Since we can 'pause' and 'start' these uploads a few items need to be considered.

### Step 1. File Selection
The first step is selecting the file you wish to upload. Once the file is selected we add it to an upload array called 'fileUploads'. Each file will also have a 'status' key that indicates where in the upload cycle it happens to be. 

Here the file should have a status of 'unauthorized'. The client should have set the status.

### Step 2. File Upload Authorization and Initialization
Once we know the file we want to upload we need to select the group it should belong to. Then we ask the server if it is ok to upload the file. The server checks the user credentials. If the credentials are good then the server creates a metadata entry in a Redis database and a file appended with 'part' in the appropriate directory. This initial file should contain 0 bytes and will be appended to as bytes are upload until the file upload is complete. 

Here the file should have a status of 'initialized'. The server should have set the status.

### Step 3. File Upload Queue
When a file is authorized we can indicate to the uploader that we wish to begin the upload of the selected file. However there can be only one file queued for upload at a time. We must loop over the files to be uploaded and if there is a previous file that has been queued we need to reset it's status to either 'paused' or 'initialized'. Whether to set the status to 'paused' or 'initialized' should be indicated by the bytes uploaded so far. If there are 0 bytes uploaded then the status should be set to 'initialized'. 

Here the file should have a status of 'queued'. The client should have set the status.

### Step 4. File Upload Start
If the file has been authorized, initialized, and queued we can begin the upload. But first we need to check if there is a current file being uploaded. If there is a file with the status of 'active' we need to send the 'pause' command to the upload worker. If there is no file currently being uploaded we can go ahead and begin or resume the file for upload.

### Step 5. File Upload Pause
If you select the pause button on a file or choose to start another file upload while there is a current upload in progress we need to pause the current upload. We set a flag in the upload worker to pause. Once the current blob upload cycle completes we check the pause flag, and if it is set to true we stop uploading blobs and tell the server that we wish to pause. All the server has to do is set the status of the appropriate file metadata to 'pause'. There should also be a file on the server with the file name and 'part' appended to the file name. That 'part' file should contain the bytes we have uploaded thus far.

Here the file should have a status of 'pause'. The server should have set the status.

Lastly once the pause cycle is complete we run the start cycle again. We check any of the upload files for a status of 'queue' (see steps 3 and 4). If a queue has been set we then begin the upload cycle for the queued file.

