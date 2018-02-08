/*
 * This class will initialize an upload to the Metis server from a web client.
 */


module.exports = function(self) {
  let timeouts = 0; // The number of times an upload has timed out.
  let maxTimeouts = 5; // The number of attepts to upload a blob.
  const INITIAL_BLOB_SIZE = Math.pow(2, 10); // in bytes

  let dispatch = (action) => self.postMessage(action);

  // report errors
  let error = (message) => dispatch(
    { type: 'WORKER_ERROR', message }
  );

  /*
   * For message passing into the worker.
   */
  let onmessage = (event) => {
    let { file, request } = event.data;

    initializeUploadSequence(file, request);
  }

  /*
   * Create an initial request to the server to initiate an upload. This lets
   * the server know what to expect and also allows the server a chance to
   * resume a previously started upload.
   */
  let initializeUploadSequence = (file, request) => {
    /*
     * Here, we are NOT attaching any blob data, but are still using a form
     * for the initialization request.
     */
    let uploadRequest = new FormData();
    for (var key in request) {
      uploadRequest.append(key,  request[key]);
    }

    /*
     * Add blob tracking information to the request. This will help us reasemble 
     * and use file upload pausing on the server. All information related to 
     * blobs is in bytes.
     */
    uploadRequest.append('current_byte_position', 0);
    uploadRequest.append('current_blob_size', 0);
    uploadRequest.append('next_blob_size', INITIAL_BLOB_SIZE);

    // Hash the next blob for security and error checking.
    let blob = file.slice(0, INITIAL_BLOB_SIZE);
    let fileReader = new FileReader();
    fileReader.onload = (progressEvent) => {
      let md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      uploadRequest.append('next_blob_hash', md5Hash);
      sendFirstPacket(uploadRequest);
    }
    fileReader.readAsArrayBuffer(blob);
  };

  let sendFirstPacket = (request) => {
    try{
      AJAX({
        url: '/upload-start',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: request,
        success: handleServerResponse,
        error: ajaxError
      });
    }
    catch({message}){
      error(message);
    }
  };

  let handleServerResponse = (response) => {
    timeouts = 0; // Reset the timeout counter;

    if (!('success' in response) || !('request' in response)) {
      error('The server response was malformed.');
      return;
    }

    if (response.request.status == 'initialized') {
      dispatch({ type: 'FILE_INITIALIZED', response });
    }
  }

  self.addEventListener('message', onmessage);
}
