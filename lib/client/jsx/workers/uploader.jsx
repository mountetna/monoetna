/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

export default (self) => {
  let dispatch = (action) => self.postMessage(action);

  let error = (message) => dispatch(
    { type: 'WORKER_ERROR', worker: 'upload', message }
  );

  // These object will be the data we need to send over the wire.
  let uploadStart = null; // The time at which the upload started in mS.
  let blobWindow = 30; // How many blob uploads to average over.
  let blobUploadTimes = [];

  let pause = false;
  let cancel = false;
  let timeouts = 0; // The number of times an upload has timed out.
  let maxTimeouts = 5; // The number of attepts to upload a blob.

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  self.addEventListener('message', ({data}) => {
    let { upload, command } = data;

    // Check that the incoming data has a valid command.
    switch (command) {
      case 'start':
        /*
         * This 'response' object usually comes as a response from
         * the server, but on a 'start' command we just fake it to
         * kick off the sequence again.
         */

        sendBlob(upload);
        break;
      case 'pause':
        // Set a pause flag.
        pause = true;
        break;
      case 'cancel':
        // Set a cancel flag.
        cancel = true;
        break;
      default:
        error('Invalid command');
        break;
    }
  });

  /*
   * Take a blob from the file and hash it.
   */
  let generateBlobHash = (blob, uploadRequest, callback) => {
    let fileReader = new FileReader();
    fileReader.onload = (event) => {
      var md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      uploadRequest.append('next_blob_hash', md5Hash);
      callback(uploadRequest);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  let sendPause = (response) => {
    var pauseData = SERIALIZE_REQUEST(response.request);
    try{
      AJAX({
        url: '/upload-pause',
        method: 'POST',
        sendType: 'serial',
        returnType: 'json',
        data: pauseData,
        success: uploader.handleServerResponse,
        error: uploader.ajaxError
      });
    }
    catch({message}) {
      error(message)
    }
  };

  let sendCancel = (response) => {
    var cancelData = SERIALIZE_REQUEST(response.request);
    try{
      AJAX({
        url: '/upload-cancel',
        method: 'POST',
        sendType: 'serial',
        returnType: 'json',
        data: cancelData,
        success: uploader.handleServerResponse,
        error: uploader.ajaxError
      });
    }
    catch(error){
      postMessage(uploader.errorObj(error.message));
    }
  }

  let sendPacket = (uploaderRequest) => {
    try{
      uploader.uploadStart = Math.floor(Date.now());
      AJAX({
        url: '/upload-blob',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: uploaderRequest,
        success: uploader.handleServerResponse,
        error: uploader.ajaxError,
        timeout: (e) => uploader.timeoutResponse(e, uploaderRequest)
      });
    }
    catch(error){
      uploader.uploadStart = null;
      postMessage(uploader.errorObj(error.message));
    }
  }

  /*
   * The request will have the squence information for the proper blob to send.
   */
  let sendBlob = (upload) => {
    let { file, current_byte_position, next_blob_size } = upload;
    /*
     * The file upload is complete, however the server should have
     * sent a 'complete' status message.
     */
    if (upload.current_byte_position >= file.size) return;

    // Recalculate the blob queue points.
    let currentBlobSize = req.next_blob_size;
    this.calcNextBlobSize(req);
    req.current_blob_size = currentBlobSize;

    // Remove old metadata from the previous packet/request
    delete req.next_blob_hash;

    // Calcuate some points on the file for chopping out blobs.
    var fromByte = req.current_byte_position;
    var toByte = req.current_byte_position + req.current_blob_size;
    var endByte = toByte + response.request.next_blob_size;

    // Set the outer limit for 'next_blob_size' and 'endByte'.
    if(endByte > fileSize){
      var cbs = req.current_blob_size;
      var cbp = req.current_byte_position;
      response.request.next_blob_size = fileSize - cbp - cbs;
      endByte = fileSize;
    }

    // Pack up a new upload packet/request
    var uploaderReq = new FormData();
    for(var key in req) uploaderReq.append(key, req[key]);

    // Append the current blob to the request.
    var blob = uploader.file.slice(fromByte, toByte);
    uploaderReq.append('blob', blob);

    // Hash the next blob.
    var nextBlob = uploader.file.slice(toByte, endByte);
    uploader.generateBlobHash(nextBlob, uploaderReq, this.sendPacket);
  };

  /* 
   * Calcuates the next blob size based upon upload speed. Also sets the last
   * calculated upload speed on the request.
   */
  let calcNextBlobSize = (request) => {
    if(uploader.uploadStart != null){
      var avgTime = uploader.averageUploadTime();
      var oldBlobSize = request.current_blob_size;
      var nextBlobSize = Math.floor(oldBlobSize * (XTR_TIME / avgTime * 0.9));

      // Set the min/max size of an upload blob.
      if(nextBlobSize > MAX_BLOB_SIZE) nextBlobSize = MAX_BLOB_SIZE;
      if(nextBlobSize < MIN_BLOB_SIZE) nextBlobSize = MIN_BLOB_SIZE;

      request.next_blob_size = nextBlobSize;

      /* 
       * Make a rough calculation of the upload speed in kilobits per second.
       * This is just for the UI.
       */
      request.uploadSpeed = (oldBlobSize * 8) / (avgTime / 1000);

      // Reset the timestamp needed for these calculations.
      uploader.uploadStart = null;
    }
  };

  let averageUploadTime = () => {
    // Get the current time it took to upload the last blob.
    var uploadTime = (Math.floor(Date.now()) - this.uploadStart);

    // Add the last upload time to an array, limit the array size with a window.
    if(this.blobUploadTimes.length >= this.blobWindow){
      this.blobUploadTimes.shift();
    }
    this.blobUploadTimes.push(uploadTime);

    // Sum the windowed upload array.
    var uploadTimeSum = this.blobUploadTimes.reduce((a, b)=>{
      return a + b;
    }, 0);

    // Average the upload time over the window.
    return uploadTimeSum / this.blobUploadTimes.length;
  }

  let handleServerResponse = (response) => {
    uploader.timeouts = 0; // Reset the timeout counter;

    var errorMessage = 'The server response was malformed.';
    if(!('success' in response) || !('request' in response)){
      postMessage(uploader.errorObj(errorMessage));
      return;
    }

    response = NORMILIZE_RESPONSE(response);

    if(!response){
      postMessage(uploader.errorObj(errorMessage));
      return;
    }

    response = TRANSFORM_RESPONSE(response);

    if(!response){
      postMessage(uploader.errorObj(errorMessage));
      return;
    }

    uploader.handleResponseRouting(response);  
  }

  let handleResponseRouting = (response) => {
    switch(response.request.status){
      case 'active':
        postMessage({ type: 'FILE_UPLOAD_ACTIVE', response });

        /*
         * Check the 'pause' and 'cancel' flags, if set then send the pause or
         * cancel message to the server.
         */
        if(uploader.pause){
          uploader.sendPause(response);
        }
        else if(uploader.cancel){
          uploader.sendCancel(response);
        }
        else{
          uploader.sendBlob(response);
        }
        break;
      case 'paused':
        /*
         * Here the server responed that it got the pause message and the upload
         * is in a paused state. We also clear the 'pause' flag on the uploader.
         */
        uploader.pause = false;
        postMessage({ type: 'FILE_UPLOAD_PAUSED', response });
        break;

      case 'cancelled':
        uploader.cancel = false;
        postMessage({ type: 'FILE_UPLOAD_CANCELLED', response });
        break;

      case 'complete':
        // Send update message back.
        uploader.file = null;
        uploader.request = null;
        postMessage({ type: 'FILE_UPLOAD_COMPLETE', response });
        break;
      default:
        // None
        break;
    }
  }

  let timeoutResponse = (error, uploadRequest) => {
    ++uploader.timeouts;
    if(uploader.timeouts == uploader.maxTimeouts){
      // Reset the uploader.
      uploader.file = null;
      uploader.request = null;
      uploader.pause = false;
      uploader.cancel = false;
      uploader.timeouts = 0;

      dispatch({ type: 'FILE_UPLOAD_TIMEOUT', response });
    }
    else{
      uploader.sendPacket(uploadRequest);
    }
  }
}
