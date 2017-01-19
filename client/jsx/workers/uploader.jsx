/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

importScripts('../libraries/utils.js');
importScripts('../libraries/spark-md5.js');

class Uploader{

  constructor(){

    // These object will be the data we need to send over the wire.
    this['file'] = null;
    this['request'] = null;

    this['uploadStart'] = null; // in milliseconds
    this['blobWindow'] = 30; // how many blob uploads to average over.
    this['blobUploadTimes'] = [];
    this['uploadSpeed'] = 0;
  }

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  onmessage(event){

    var message = event['data'];
    uploader['file'] = message['file'];
    uploader['request'] = message['request'];

    // Check that the incoming data has a valid command.
    if(!'command' in message){

      postMessage({ message: 'No command.'});
      return;
    }

    var command = message['command'];
    switch(command){

      case 'initialize':

        uploader.initializeUploadSequence();
        break;
      case 'start':

        /* 
         * This 'response' object usually comes as a response from the server, 
         * but on a 'start' command we just fake it to kick of the sequence
         * again.
         */

        uploader.snakeCaseIt(uploader['request']);

        var response = {

          'request': uploader['request'],
          'byte_count':  uploader['request']['current_byte_position']
        };

        uploader.sendBlob(response);
        break;
      case 'pause':

        // Pause the upload.
        break;
      case 'query':

        // Query the server for the current upload status.
        break;
      default:

        break;
    }
  }

  /*
   * Create an initial request to the server to initiate an upload. This lets
   * the server know what to expect and also allows the server a chance to
   * resume a previously started upload.
   */
  initializeUploadSequence(){

    /*
     * Here, we are NOT attaching any blob data, but are still using a form
     * for the initialization request. Also keep an eye on that 'SNAKE_CASE_IT'
     * function. The Thin/Ruby server wants it's vars in snake case.
     */
    var initUpReq = new FormData();
    for(var key in  uploader['request']){

      initUpReq.append(SNAKE_CASE_IT(key),  uploader['request'][key]);
    }

    /*
     * Add blob tracking information to the request. This will help us reasemble 
     * and use file upload pausing on the server. All information related to 
     * blobs is in bytes.
     */
    initUpReq.append('current_byte_position', 0);
    initUpReq.append('current_blob_size', 0);
    initUpReq.append('next_blob_size', BLOB_SIZE);

    // Hash the next blob for security and error checking.
    var blob = uploader['file'].slice(0, BLOB_SIZE);
    uploader.generateBlobHash(blob, initUpReq, this.sendFirstPacket);
  }

  /*
   * Take a blob from the file and hash it.
   */
  generateBlobHash(blob, uploadRequest, callback){

    var fileReader = new FileReader();

    fileReader.onload = function(progressEvent){

      var md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      uploadRequest.append('next_blob_hash', md5Hash);

      callback(uploadRequest);
    }

    fileReader.readAsArrayBuffer(blob);
  }

  sendFirstPacket(initialUploadRequet){

    try{

      AJAX({

        'url': '/upload-start',
        'method': 'POST',
        'sendType': 'file',
        'returnType': 'json',
        'data': initialUploadRequet,
        'success': uploader['handleServerResponse'],
        'error': uploader['ajaxError']
      });
    }
    catch(error){

      uploader['uploadStart'] = null;
      postMessage({ 'type': 'error', 'message': error['message'] });
    }
  }

  /*
   * The request will have the squence information for the proper blob to send.
   */
  sendBlob(response){

    uploader.generateUploaderRequest(response);
  }

  generateUploaderRequest(response){

    var fileSize = uploader['file'].size
    if(response['byte_count'] >= fileSize){

      /*
       * The file upload is complete, however the server should have sent a 
       * 'complete' status message.
       */
      return;
    }

    var request = response['request'];

    // Set the blob cue points.
    var fromByte = parseInt(response['byte_count']);
    var toByte = fromByte + BLOB_SIZE;
    toByte = (toByte > fileSize) ? fileSize : toByte;

    /*
     * Set the the metadata for the next blob in sequence. This part gets used 
     * by the server for data integrety and a wee bit of security.
     */
    request['current_byte_position'] = fromByte;
    request['current_blob_size'] = toByte - fromByte;

    // Remove old metadata from the previous packet/request
    delete request['next_blob_hash'];
    delete request['next_blob_size'];

    // Pack up a new upload packet/request
    var uploaderRequest = new FormData();
    for(var key in request){

      uploaderRequest.append(key, request[key]);
    }

    // Append the current blob to the request.
    var blob = uploader['file'].slice(fromByte, toByte);
    uploaderRequest.append('blob', blob);

    // Hash and profile the next blob for security and error checking.
    var endByte = toByte + NEXT_BLOB_SIZE;
    endByte = (endByte > fileSize) ? fileSize : endByte;
    uploaderRequest.append('next_blob_size', (endByte - toByte));
    var nextBlob = uploader['file'].slice(toByte, endByte);
    uploader.generateBlobHash(nextBlob, uploaderRequest, this.sendPacket);
  }

  sendPacket(uploaderRequest){

    try{

      uploader['uploadStart'] = Math.floor(Date.now());

      AJAX({

        'url': '/upload-blob',
        'method': 'POST',
        'sendType': 'file',
        'returnType': 'json',
        'data': uploaderRequest,
        'success': uploader['handleServerResponse'],
        'error': uploader['ajaxError']
      });
    }
    catch(error){

      uploader['uploadStart'] = null;
      postMessage({ type: 'error', message: error['message'] });
    }
  }

  // Set blob size based upon upload time.
  setBlobSize(){

    if(uploader['uploadStart'] != null){

      var avgTime = uploader.averageUploadTime();
      var olbBlobSize = BLOB_SIZE;
      BLOB_SIZE = NEXT_BLOB_SIZE;
      NEXT_BLOB_SIZE = Math.floor(olbBlobSize * (TRANSFER_TIME / avgTime * 0.9))

      uploader['uploadStart'] = null;
      if(NEXT_BLOB_SIZE < MIN_BLOB_SIZE) NEXT_BLOB_SIZE = MIN_BLOB_SIZE;

      // Make a rough calculation of the upload speed in kilobits per second.
      this['uploadSpeed'] = (BLOB_SIZE * 8) / (avgTime/1000);
    }
  }

  averageUploadTime(){

    // Get the current time it took to upload the last blob.
    var uploadTime = (Math.floor(Date.now()) - uploader['uploadStart']);

    // Add the last upload time to an array, limit the array size with a window.
    if(this['blobUploadTimes']['length'] >= this['blobWindow']){

      this['blobUploadTimes'].shift();
    }
    this['blobUploadTimes'].push(uploadTime);

    // Sum the windowed upload array.
    var uploadTimeSum = this['blobUploadTimes'].reduce((a, b)=>{

      return a + b;
    }, 0);

    // Average the upload time over the window.
    return uploadTimeSum / this['blobUploadTimes']['length'];
  }

  handleServerResponse(response){

    // Set blob size based upon upload time.
    uploader.setBlobSize();

    var errorMessage = 'The server response was malformed.';
    if(!('success' in response) || !('status' in response)){

      postMessage({ 'type': 'error', message: errorMessage });
      return;
    }

    response = NORMILIZE_RESPONSE(response);

    if(!response){

      postMessage({ type: 'error', message: errorMessage });
      return;
    }

    response = TRANSFORM_RESPONSE(response);

    if(!response){

      postMessage({ type: 'error', message: errorMessage });
      return;
    }

    uploader.handleResponseRouting(response);  
  }

  handleResponseRouting(response){

    switch(response['status']){

      case 'initialized':

        postMessage({ type: 'initialized', response: response });
        break;
      case 'active':

        response['request']['uploadSpeed'] = this['uploadSpeed'];
        postMessage({ type: 'active', response: response });
        uploader.sendBlob(response);
        break;
      case 'paused':

        break;
      case 'complete':

        // Send update message back.
        this['file'] = null;
        this['request'] = null;
        postMessage({ type: 'complete', response: response });
        break;
      case 'stopped':

        break;
      case 'failed':

        break;
      default:

        // None
        break;
    }
  }

  snakeCaseIt(object){

    for(var key in object){

      object[SNAKE_CASE_IT(key)] = object[key];
      if(key != key.toLowerCase()) delete object[key];
    }
  }

  /*
   * Handle the AJAX errors.
   */
  ajaxError(xhr, ajaxOpt, thrownError){

    postMessage({ xhr: xhr, ajaxOptions: ajaxOpt, thrownError: thrownError });
  }
}

// Initilize the class.
var uploader = new Uploader();

// Expose the worker specific messaging function.
var onmessage = uploader.onmessage;