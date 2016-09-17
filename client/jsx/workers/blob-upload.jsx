/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

importScripts('../utils.js');

class BlobUpload{

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  onmessage(event){

    var message = event['data'];

    var request = message['data']['request'];
    var signature = message['data']['signature'];
    
    blobUpload['file'] = message['data']['uploadFile'];

    // Check that the incoming data has a valid command.
    if(!'command' in message){

      postMessage({ message: 'No command.'});
      return;
    }

    var command = message['command'];
    switch(command){

      case 'start':

        
        blobUpload.initializeUploadSequence(request, signature);
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

    var initUpReq = blobUpload.generateInitialUploadRequest(req, sig);

    try{

      AJAX({

        url: '/upload-init',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: initUpReq,
        success: this.handleServerResponse,
        error: this.ajaxError
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

    var blobUpReq = blobUpload.generateBlobUploadRequest(response);

    try{

      AJAX({

        url: '/upload-blob',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: blobUpReq,
        success: this.handleServerResponse,
        error: this.ajaxError
      });
    }
    catch(error){

      postMessage({ 'type': 'error', message: error['message'] });
    }
  }

  generateInitialUploadRequest(request, signature){

    // Append the authorization signature/hash from Magma to the request
    request['signature'] = signature;

    /*
     * Add blob tracking information to the request. This will help us reasemble 
     * and use file upload pausing on the server. All information related to 
     * blobs is in bytes.
     */
    request['current_byte_position'] = 0;
    request['current_blob_size'] = 0;

    request['next_blob_size'] = BLOB_SIZE;
    request['next_blob_hash'] = blobUpload.generateBlobHash(0, BLOB_SIZE);

    /*
     * Here, we are NOT attaching any blob data, but are still using the form
     * for the initialization request.
     */
    var initialUploadRequest = new FormData();
    for(var key in request){

      initialUploadRequest.append(key, request[key]);
    }

    return initialUploadRequest;
  }

  generateBlobUploadRequest(response){

    var request = response['request'];

    // Set the blob cue points.
    var fromByte = parseInt(response['byte_count']);
    var toByte = parseInt(fromByte) + BLOB_SIZE;
    
    request['current_byte_position'] = fromByte;
    request['current_blob_size'] = BLOB_SIZE;

    /*
     * Set the the metadata for the next blob in sequence. This bit gets used by
     * the server for data integrety and a wee bit of security.
     */
    var endByte = toByte + BLOB_SIZE;
    request['next_blob_hash'] = blobUpload.generateBlobHash(toByte, endByte);
    request['next_blob_size'] = BLOB_SIZE;

    // Slice out a blob from our upload file.
    var blob = blobUpload['file'].slice(fromByte, toByte);

    // Pack all of our data up into a FormData object
    var blobUploadRequest = new FormData();
    blobUploadRequest.append('blob', blob);

    for(var key in request){

      blobUploadRequest.append(key, request[key]);
    }

    return blobUploadRequest;
  }

  /*
   * Take a blob from the file and hash it.
   */
  generateBlobHash(fromByte, toByte){

    /*
     * We need to insert a proper hashing algo here, but for now lets just 
     * send back a dummy hash.
     */

    // use the singleton variable blobUpload['file'];
    return 'blahf8cfc63531b6ed753bd536f6e12d578c';
  }

  handleServerResponse(response){

    var errorMessage = 'The server did not give a valid response.'

    response = VERIFY_AND_TRANSFORM(response);
    if(!response){
  
      postMessage({ 'type': 'error', message: errorMessage });
    }

    switch(response['status']){

      case 'initialized':
      case 'active':

        blobUpload.sendBlob(response);
        break;
      case 'paused':

        break;
      case 'complete':

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
var blobUpload = new BlobUpload();

// Expose the worker specific messaging function.
var onmessage = blobUpload.onmessage;