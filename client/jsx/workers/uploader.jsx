/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

importScripts('../lib/utils.js');
importScripts('../lib/spark-md5.js');

class Uploader{

  constructor(){

    // This object will be the data we need to send over the wire.
    this['file'] = null; 
  }

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  onmessage(event){

    var message = event['data'];

    var request = message['data']['request'];
    var signature = message['data']['signature'];
    
    uploader['file'] = message['data']['uploadFile'];

    // Check that the incoming data has a valid command.
    if(!'command' in message){

      postMessage({ message: 'No command.'});
      return;
    }

    var command = message['command'];
    switch(command){

      case 'start':
        
        uploader.initializeUploadSequence(request, signature);
        break;

      case 'pause':

        // Pause the upload.
        break;
      case 'stop':

        // Stop the upload.
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
  initializeUploadSequence(req, sig){

    //this.initUpReq = uploader.generateInitialUploadRequest(req, sig);
    uploader.generateInitialUploadRequest(req, sig);
  }

  generateInitialUploadRequest(request, signature){

    // Append the authorization signature/hash from Magma to the request
    request['signature'] = signature;

    /*
     * Here, we are NOT attaching any blob data, but are still using the form
     * for the initialization request.
     */
    var initUpReq = new FormData();
    for(var key in request){

      initUpReq.append(key, request[key]);
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

        url: '/upload-start',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: initialUploadRequet,
        success: uploader.handleServerResponse,
        error: uploader.ajaxError
      });
    }
    catch(error){

      postMessage({ 'type': 'error', message: error['message'] });
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
    var endByte = toByte + BLOB_SIZE;
    endByte = (endByte > fileSize) ? fileSize : endByte;
    uploaderRequest.append('next_blob_size', (endByte - toByte));
    var nextBlob = uploader['file'].slice(toByte, endByte);
    uploader.generateBlobHash(nextBlob, uploaderRequest, this.sendPacket);
  }

  sendPacket(uploaderRequest){

    try{

      AJAX({

        url: '/upload-blob',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: uploaderRequest,
        success: uploader.handleServerResponse,
        error: uploader.ajaxError
      });
    }
    catch(error){

      postMessage({ 'type': 'error', message: error['message'] });
    }
  }

  handleServerResponse(response){

    /*
    if(response === false){
    
      var errorMessage = 'The server did not give a valid response.';
      postMessage({ 'type': 'error', message: errorMessage });
    }
    */

    switch(response['status']){

      case 'initialized':
      case 'active':

        response = VERIFY_AND_TRANSFORM(response);
        uploader.sendBlob(response);
        break;
      case 'paused':

        break;
      case 'complete':

        console.log('Complete.');
        console.log(response);
        break;
      case 'stopped':

        break;
      case 'failed':

        break;
      default:

        //none
        break;
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